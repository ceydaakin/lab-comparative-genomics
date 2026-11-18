"""Collect finished TYGS jobs (tygs.dsmz.de; web form with e-mail, no submission API) by job GUID and summarise dDDH.
d4 (GGDC formula 2) is the TYGS-recommended value for draft genomes; >= 70% = same species. A d4 >= 70% to another genus
with a low d0 (whole-genome formula) flags a few co-binned foreign contigs in a MAG, not a species question."""
import argparse, io, os, urllib.error, urllib.request
import numpy as np, pandas as pd
from parsers import html_text

ap = argparse.ArgumentParser()
ap.add_argument("--jobs", nargs="*", default=[]); ap.add_argument("--genomes", required=True)
ap.add_argument("--metadata", required=True); ap.add_argument("--outdir", required=True)
A = ap.parse_args()
os.makedirs(A.outdir, exist_ok=True)
G = [l.strip() for l in open(A.genomes) if l.strip()]
meta = pd.read_csv(A.metadata, sep="\t").set_index("accession")
BASE = "https://tygs.dsmz.de/user_results"



def fetch(url, guid):
    try:
        return urllib.request.urlopen(url, timeout=300).read()
    except urllib.error.HTTPError as e:
        raise SystemExit(f"TYGS job {guid}: HTTP {e.code} for {url} (wrong GUID, expired job, or TYGS down)")


frames = []
for guid in A.jobs:
    table = pd.read_csv(io.BytesIO(fetch(f"{BASE}/show/digital_ddh_tsv/{guid}", guid)), sep="\t")
    if "comp_type" not in table.columns or table.empty:
        raise SystemExit(f"TYGS job {guid} has no dDDH rows (job not finished, or the GUID is wrong)")
    frames.append(table.assign(job=guid))
    with open(f"{A.outdir}/TYGS_report_{guid[:8]}.pdf", "wb") as f:
        f.write(fetch(f"{BASE}/show.pdf?guid={guid}", guid))

cols = ["genome", "strain", "type_strain", "digital_ddh_d4", "confidence_interval_d4", "digital_ddh_d0", "species_call", "foreign_type_strain_d4>=70"]
if not frames:
    pd.DataFrame(columns=cols).to_csv(f"{A.outdir}/tygs_species_identification.tsv", sep="\t", index=False)
    open(f"{A.outdir}/coverage.txt", "w").write(f"0/{len(G)} genomes in TYGS jobs (none configured)\n")
    raise SystemExit(0)

d = pd.concat(frames)
d["query"], d["subject"] = d.query_genome.map(html_text), d.subject_genome.map(html_text)
d = d.drop(columns=["query_genome", "subject_genome"])
d.to_csv(f"{A.outdir}/tygs_digital_ddh.tsv", sep="\t", index=False)
ut = d[(d.comp_type == "U vs. T") & (d["query"].isin(G) | d["subject"].isin(G))].copy()   # a job may hold extra genomes
ut["genome"] = [q if q in G else s for q, s in zip(ut["query"], ut["subject"])]
ut["type_strain"] = [s if q in G else q for q, s in zip(ut["query"], ut["subject"])]
best = ut.sort_values("digital_ddh_d4", ascending=False).groupby("genome").head(1).copy()
best["strain"] = best.genome.map(meta.strain)
best["species_call"] = np.where(best.digital_ddh_d4 >= 70, "same species as type strain", "potential novel species")
genus = best.type_strain.str.split().str[0]
foreign = ut[(ut.type_strain.str.split().str[0] != ut.genome.map(dict(zip(best.genome, genus)))) & (ut.digital_ddh_d4 >= 70)]
best["foreign_type_strain_d4>=70"] = best.genome.map(foreign.groupby("genome").type_strain.agg("; ".join)).fillna("")
best[cols].sort_values("digital_ddh_d4").to_csv(f"{A.outdir}/tygs_species_identification.tsv", sep="\t", index=False)
uu = d[d.comp_type == "U vs. U"]
m = pd.DataFrame(100.0, index=G, columns=G)
for q, s, v in zip(uu["query"], uu["subject"], uu.digital_ddh_d4):
    if q in m.index and s in m.columns: m.loc[q, s] = m.loc[s, q] = v
m.to_csv(f"{A.outdir}/tygs_dddh_d4_matrix.tsv", sep="\t")
covered = set(best.genome) & set(G)
open(f"{A.outdir}/coverage.txt", "w").write(f"{len(covered)}/{len(G)} genomes in TYGS jobs; missing: {sorted(set(G) - covered)}\n")
print(f"TYGS: {len(covered)}/{len(G)} genomes; species calls {best.species_call.value_counts().to_dict()}")
