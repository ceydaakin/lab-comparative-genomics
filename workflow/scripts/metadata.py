"""Genome metadata from the NCBI Datasets assembly report (JSON Lines) plus a niche label per genome.
Niche = first config keyword found in isolation_source / geo / strain (case-insensitive), else the default label."""
import argparse, csv, json

p = argparse.ArgumentParser()
p.add_argument("--report", required=True)
p.add_argument("--niches", required=True, help="JSON list of [keyword, label] pairs")
p.add_argument("--default-niche", default="Other")
p.add_argument("--out", required=True)
a = p.parse_args()
niches = json.loads(a.niches)

rows = []
for line in open(a.report):
    d = json.loads(line)
    ai = d.get("assemblyInfo", {}); bs = ai.get("biosample", {})
    attrs = {x["name"]: x.get("value", "") for x in bs.get("attributes", [])}
    inf = d.get("organism", {}).get("infraspecificNames", {})
    st, ck = d.get("assemblyStats", {}), d.get("checkmInfo", {})
    source = bs.get("isolationSource") or attrs.get("isolation_source", "")
    strain = inf.get("strain") or inf.get("isolate") or bs.get("strain") or attrs.get("strain", "") or d["accession"]
    text = " ".join([source, bs.get("geoLocName", ""), strain]).lower()
    niche = next((label for kw, label in niches if kw.lower() in text), a.default_niche)
    rows.append({"accession": d["accession"], "strain": strain, "niche": niche, "source": source,
                 "geo": bs.get("geoLocName") or attrs.get("geo_loc_name", ""),
                 "date": bs.get("collectionDate") or attrs.get("collection_date", ""),
                 "level": ai.get("assemblyLevel"), "size": st.get("totalSequenceLength"),
                 "contigs": st.get("numberOfContigs"), "n50": st.get("contigN50"), "gc": st.get("gcPercent"),
                 "checkm_completeness": ck.get("completeness"), "checkm_contamination": ck.get("contamination"),
                 "bioproject": ai.get("bioprojectAccession", "")})
with open(a.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
