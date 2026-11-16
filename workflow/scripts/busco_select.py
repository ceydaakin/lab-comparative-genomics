"""Collect BUSCO short-summary JSONs, apply the completeness cut-off and write the analysis genome set.
Fails loudly if the focal genome does not pass, or if fewer than 3 genomes remain."""
import argparse, json, os
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--summaries", nargs="+", required=True)
p.add_argument("--metadata", required=True)
p.add_argument("--min-complete", type=float, required=True)
p.add_argument("--exclude", nargs="*", default=[])
p.add_argument("--focal", required=True)
p.add_argument("--table", required=True)
p.add_argument("--set", required=True)
a = p.parse_args()

rows = []
for f in a.summaries:
    acc = os.path.basename(f).removesuffix(".json")   # results/qc/busco/<accession>.json
    r = json.load(open(f))["results"]
    rows.append({"accession": acc, "busco": r["one_line_summary"], "complete": float(r["Complete percentage"]),
                 "duplicated": float(r["Multi copy percentage"]), "fragmented": float(r["Fragmented percentage"]),
                 "missing": float(r["Missing percentage"])})
t = pd.DataFrame(rows).merge(pd.read_csv(a.metadata, sep="\t"), on="accession", how="left")
t["excluded_by_config"] = t.accession.isin(a.exclude)
t["included"] = (t.complete >= a.min_complete) & ~t.excluded_by_config
t["reason"] = ["listed in config qc.exclude" if x else (f"BUSCO {c}% < {a.min_complete}%" if not i else "")
               for x, c, i in zip(t.excluded_by_config, t.complete, t.included)]
t.to_csv(a.table, sep="\t", index=False)
keep = t[t.included].accession.tolist()
if a.focal not in keep:
    raise SystemExit(f"focal genome {a.focal} failed QC ({t.set_index('accession').reason.get(a.focal, 'not in input')})")
if len(keep) < 3:
    raise SystemExit(f"only {len(keep)} genomes passed QC; at least 3 are needed")
open(a.set, "w").write("\n".join(keep) + "\n")
print(f"{len(keep)}/{len(t)} genomes pass BUSCO >= {a.min_complete}%")
