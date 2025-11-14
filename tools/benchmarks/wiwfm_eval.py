#!/usr/bin/env python3
import argparse, csv, json, time
import httpx

# Simple, conservative Yes/No policy
# - YES if tier == I and confidence >= 0.75 and no resistance flags
# - NO if resistance:* present OR confidence < 0.4
# - ELSE MAYBE (research)

def has_resistance(rationale):
    return any(str(x).startswith("resistance:") for x in (rationale or []))


def decide(pred):
    tier = pred.get("tier")
    conf = float(pred.get("confidence") or 0.0)
    rat = pred.get("rationale") or []
    if has_resistance(rat) or conf < 0.4:
        return "NO"
    if tier == "I" and conf >= 0.75 and not has_resistance(rat):
        return "YES"
    return "MAYBE"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--api_base", default="http://127.0.0.1:8000")
    ap.add_argument("--mode", default="chemo", choices=["chemo","radonc","synthetic"])
    ap.add_argument("--timeout", type=float, default=60.0)
    ap.add_argument("--retries", type=int, default=2)
    ap.add_argument("--max_rows", type=int, default=0, help="limit number of rows to process; 0 means all")
    ap.add_argument("--sleep", type=float, default=0.15, help="per-row sleep throttle seconds")
    args = ap.parse_args()

    rows=[]
    with open(args.csv, newline="", encoding="utf-8") as f:
        rows_all = list(csv.DictReader(f))
        rows = rows_all[: args.max_rows] if args.max_rows and args.max_rows > 0 else rows_all

    out=[]
    with httpx.Client(timeout=args.timeout) as client:
        for r in rows:
            payload={"mutations":[], "api_base": args.api_base}
            disease=r.get("disease")
            therapy=r.get("therapy") or r.get("drug_or_class")
            gene=r.get("gene"); hg=r.get("hgvs_p")
            chrom=r.get("chrom"); pos=r.get("pos"); ref=r.get("ref"); alt=r.get("alt"); build=r.get("build")
            t_lc = (therapy or "").strip().lower()
            options = {"adaptive": True, "ensemble": True, "disable_literature": True}
            if args.mode=="chemo":
                # Route platinum/PARP to synthetic automatically
                if t_lc in {"platinum", "parp inhibitor", "parp-inhibitor"}:
                    endpoint=f"{args.api_base}/api/guidance/synthetic_lethality"; payload.update({"disease":disease, "options": options})
                else:
                    endpoint=f"{args.api_base}/api/guidance/chemo"; payload.update({"disease":disease,"drug_or_class":therapy, "options": options})
            elif args.mode=="radonc":
                endpoint=f"{args.api_base}/api/guidance/radonc"; payload.update({"disease":disease, "options": options})
            else:
                endpoint=f"{args.api_base}/api/guidance/synthetic_lethality"; payload.update({"disease":disease, "options": options})
            mut={"gene":gene,"hgvs_p":hg}
            if chrom and pos and ref and alt:
                mut.update({"chrom":str(chrom),"pos":int(pos),"ref":str(ref),"alt":str(alt)})
            if build:
                mut["build"]=build
            payload["mutations"]=[mut]
            pred={}
            last_err=None
            for attempt in range(max(1, args.retries)):
                try:
                    resp=client.post(endpoint,json=payload)
                    if resp.status_code<400:
                        pred=resp.json() or {}
                        break
                    else:
                        last_err=resp.text
                except Exception as e:
                    last_err=str(e)
                time.sleep(0.4 * (attempt+1))
            if not pred:
                pred={"error": last_err or "timed out"}
            decision = decide(pred) if "error" not in pred else "ERROR"
            out.append({"input": r, "prediction": pred, "decision": decision})
            time.sleep(max(0.0, args.sleep))
    print(json.dumps(out, indent=2))

if __name__=="__main__":
    main()

            last_err=None
            for attempt in range(max(1, args.retries)):
                try:
                    resp=client.post(endpoint,json=payload)
                    if resp.status_code<400:
                        pred=resp.json() or {}
                        break
                    else:
                        last_err=resp.text
                except Exception as e:
                    last_err=str(e)
                time.sleep(0.4 * (attempt+1))
            if not pred:
                pred={"error": last_err or "timed out"}
            decision = decide(pred) if "error" not in pred else "ERROR"
            out.append({"input": r, "prediction": pred, "decision": decision})
            time.sleep(max(0.0, args.sleep))
    print(json.dumps(out, indent=2))

if __name__=="__main__":
    main()
