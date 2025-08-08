
import modal
import pytest

# A sample protein sequence for testing. This is a real, short protein segment.
# It's part of the Human Serum Albumin.
SAMPLE_PROTEIN = "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL"

@pytest.mark.integration
def test_command_center_esm_sieve():
    """
    Tests the isolated ESM Sieve endpoint on the CommandCenter.
    """
    # FIX: Use the correct, deployed app name for lookup
    CommandCenter = modal.Cls.lookup("command-center-v8-override-fix", "CommandCenter")
    cc = CommandCenter()
    
    print("\n--- 🔬 Initiating isolated ESM Sieve test... ---")
    
    score = cc.test_esm_sieve.remote(SAMPLE_PROTEIN)
    
    assert isinstance(score, float), f"Expected a float score, but got {type(score)}"
    assert score != 0.0, "Expected a non-zero score." # It's highly unlikely to be exactly 0.0
    
    print(f"--- ✅ SUCCESS: Received ESM Score: {score:.4f} ---")

if __name__ == "__main__":
    test_command_center_esm_sieve() 