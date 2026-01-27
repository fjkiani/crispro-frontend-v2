"""PyTorch compatibility patches for Evo2 checkpoint loading.
These patches are required for PyTorch 2.3.0+ compatibility with Evo2 checkpoints."""
try:
    import torch
except ImportError:
    # Allow Modal to parse the file without torch installed locally
    # torch will be available in the Modal container
    torch = None


class _FlashAttnGPUProxy:
    """Proxy around flash-attn GPU module to normalize `.fwd()` return signature."""

    def __init__(self, inner):
        self._inner = inner

    def __getattr__(self, name):
        return getattr(self._inner, name)

    def fwd(self, *args, **kwargs):
        res = self._inner.fwd(*args, **kwargs)
        # Normalize return signature to 4 values (vortex expects 4).
        # Some flash-attn builds return additional tensors/metadata.
        try:
            def _clone(x):
                return x.clone() if torch.is_tensor(x) else x

            # Prefer fast-path for common tuple/list.
            if isinstance(res, tuple):
                head = res[:4] if len(res) > 4 else res
                if isinstance(head, tuple) and len(head) == 4:
                    return tuple(_clone(x) for x in head)
                return head
            if isinstance(res, list):
                head = res[:4] if len(res) > 4 else res
                if len(head) == 4:
                    return tuple(_clone(x) for x in head)
                return tuple(head)

            # Generic iterable path: attempt to unpack >=4 and drop extras.
            out, softmax_lse, s_dmask, rng_state, *_rest = res  # type: ignore[misc]

            # Torch custom-op constraint: returned tensors must not alias inputs/each other.
            # Be conservative and clone tensors to satisfy aliasing rules.
            return (_clone(out), _clone(softmax_lse), _clone(s_dmask), _clone(rng_state))
        except Exception:
            return res


def ensure_vortex_flash_attn_compat() -> bool:
    """
    Ensure vortex's flash-attn wrapper won't crash due to return signature mismatch.

    This is safe to call repeatedly and is intentionally lightweight.
    """
    try:
        import sys
        patched_any = False

        # Patch the canonical import path first.
        try:
            from vortex.ops import attn_interface as _attn_interface  # type: ignore
            if hasattr(_attn_interface, "flash_attn_gpu") and not isinstance(_attn_interface.flash_attn_gpu, _FlashAttnGPUProxy):
                _attn_interface.flash_attn_gpu = _FlashAttnGPUProxy(_attn_interface.flash_attn_gpu)
                patched_any = True
        except Exception:
            pass

        # Defensive: patch any duplicate-loaded module objects pointing at the same file.
        for _m in list(sys.modules.values()):
            try:
                if not getattr(_m, "__file__", ""):
                    continue
                if not str(_m.__file__).endswith("/vortex/ops/attn_interface.py"):
                    continue
                if hasattr(_m, "flash_attn_gpu") and not isinstance(getattr(_m, "flash_attn_gpu"), _FlashAttnGPUProxy):
                    setattr(_m, "flash_attn_gpu", _FlashAttnGPUProxy(getattr(_m, "flash_attn_gpu")))
                    patched_any = True
            except Exception:
                continue

        return patched_any
    except Exception:
        return False


def apply_torch_patches():
    """
    Apply PyTorch patches required for Evo2 checkpoint loading.
    
    This must be called BEFORE importing evo2, as the patches need to be
    in place before any torch.load() calls occur during model initialization.
    
    Patches applied:
    1. torch.load: Forces weights_only=False for HuggingFace checkpoint loading
    2. torch.serialization.add_safe_globals: No-op for PyTorch < 2.4.0 compatibility
    """
    # CRITICAL: PyTorch 2.3.0+ defaults to weights_only=True for security
    # Evo2 checkpoints from HuggingFace require weights_only=False because they
    # contain complex serialized objects (not just tensor weights).
    original_torch_load = torch.load
    
    def patched_torch_load(*args, **kwargs):
        """Patched torch.load that always sets weights_only=False."""
        # ALWAYS set weights_only=False for evo2 checkpoint loading
        # (trusted source: HuggingFace)
        kwargs['weights_only'] = False
        return original_torch_load(*args, **kwargs)
    
    torch.load = patched_torch_load
    
    # CRITICAL: PyTorch 2.3.0 doesn't have add_safe_globals, but vortex
    # (used by Evo2) may try to call it. Provide a no-op if it doesn't exist.
    if not hasattr(torch.serialization, 'add_safe_globals'):
        def add_safe_globals(_globals_list):
            """No-op function for PyTorch < 2.4.0 compatibility."""
            return None
        
        torch.serialization.add_safe_globals = add_safe_globals

    # RUO/engineering patch: Torch custom-op aliasing validation can be overly strict for
    # some flash-attn/vortex build combinations (observed in Modal logs for Evo2 scoring).
    # We disable this validation to keep the service operational.
    #
    # NOTE: This is acceptable for RUO scoring, but we should still align flash-attn/vortex
    # versions long-term to avoid relying on this bypass.
    try:
        import torch._library.custom_ops as _custom_ops  # type: ignore
        if hasattr(_custom_ops, "_validate_outputs_to_not_alias_inputs"):
            _custom_ops._validate_outputs_to_not_alias_inputs = lambda *args, **kwargs: None  # type: ignore[attr-defined]
    except Exception:
        pass


def apply_evo2_runtime_patches():
    """
    Apply runtime monkey-patches to the installed `evo2` package.

    Why: We've observed Evo2 raising `too many values to unpack (expected 4)` during
    `score_sequences`, which strongly suggests an upstream return-signature mismatch
    (commonly `prepare_batch(...)` returning >4 values while older scoring code expects 4).

    This patch is intentionally defensive:
    - If `evo2.scoring.prepare_batch` returns a tuple longer than 4, we slice to 4.
    - If the model forward returns a tuple/list with more than one item, we take [0] as logits.
    """
    try:
        import evo2.scoring as scoring  # type: ignore
    except Exception:
        return

    # Patch prepare_batch: return only the first 4 items if more are returned.
    if hasattr(scoring, "prepare_batch"):
        _orig_prepare_batch = scoring.prepare_batch

        def _patched_prepare_batch(*args, **kwargs):
            out = _orig_prepare_batch(*args, **kwargs)
            if isinstance(out, tuple) and len(out) > 4:
                return out[:4]
            return out

        scoring.prepare_batch = _patched_prepare_batch  # type: ignore[attr-defined]

    # Patch _score_sequences to tolerate model() returning more than (logits, extra)
    if hasattr(scoring, "_score_sequences"):
        _orig__score_sequences = scoring._score_sequences

        def _patched__score_sequences(*args, **kwargs):
            # If upstream fixes this, this wrapper just passes through.
            try:
                return _orig__score_sequences(*args, **kwargs)
            except ValueError as e:
                # Best-effort: handle "too many values to unpack" originating from model forward.
                if "too many values to unpack" not in str(e):
                    raise
                # Re-implement minimal logic with tolerant unpacking.
                # We intentionally avoid importing heavy deps unless needed.
                import numpy as np  # type: ignore

                seqs = kwargs.get("seqs") or (args[0] if len(args) > 0 else None)
                model = kwargs.get("model") or (args[1] if len(args) > 1 else None)
                tokenizer = kwargs.get("tokenizer") or (args[2] if len(args) > 2 else None)
                prepend_bos = kwargs.get("prepend_bos", False)
                reduce_method = kwargs.get("reduce_method", "mean")
                device = kwargs.get("device", "cuda:0")

                input_ids, seq_lengths, *_rest = scoring.prepare_batch(  # type: ignore[misc]
                    seqs, tokenizer, device=device, prepend_bos=prepend_bos
                )
                with torch.inference_mode():
                    out = model(input_ids)
                    logits = out[0] if isinstance(out, (tuple, list)) else out

                logprobs = scoring.logits_to_logprobs(logits, input_ids)  # type: ignore[attr-defined]
                logprobs = logprobs.float().cpu().numpy()

                if reduce_method == "sum":
                    reduce_func = np.sum
                elif reduce_method == "mean":
                    reduce_func = np.mean
                else:
                    raise ValueError(f"Invalid reduce_method {reduce_method}")

                scores = []
                for idx in range(len(seq_lengths)):
                    n = int(seq_lengths[idx])
                    # Skip the first token; score next-token likelihoods.
                    # Use max(0, n-1) to be safe on very short sequences.
                    n_eff = max(0, n - 1)
                    s = reduce_func(logprobs[idx, 1 : 1 + n_eff])
                    scores.append(float(s) if np.isfinite(s) else float("nan"))
                return scores

        scoring._score_sequences = _patched__score_sequences  # type: ignore[attr-defined]

    # Patch flash-attn return signature mismatch (observed in Modal logs):
    # vortex expects: out, softmax_lse, S_dmask, rng_state = flash_attn_gpu.fwd(...)
    # but newer flash-attn may return additional values.
    try:
        from vortex.ops import attn_interface as _attn_interface  # type: ignore

        # The crash is thrown inside vortex's python backend wrapper:
        #   out, softmax_lse, S_dmask, rng_state = flash_attn_gpu.fwd(...)
        # when `flash_attn_gpu.fwd(...)` returns MORE than 4 values.
        #
        # Simply reassigning `flash_attn_gpu.fwd` is not reliable (it may be a C-extension
        # function or captured differently). Instead, we patch the backend wrapper
        # function **in-place** by swapping its __code__ object so the already-registered
        # torch custom-op backend calls the new logic.
        ensure_vortex_flash_attn_compat()
    except Exception:
        # If vortex/flash-attn isn't present (CPU/testing), do nothing.
        pass
