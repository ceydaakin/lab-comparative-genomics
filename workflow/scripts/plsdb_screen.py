"""Parse mash screen -w output (one file per genome) against PLSDB and annotate hits with plasmid metadata.
A plasmid is called 'present' when at least --present-fraction of its sketch hashes are shared."""
import argparse, os
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--screens", nargs="+", required=True); p.add_argument("--nuccore", required=True); p.add_argument("--taxonomy", required=True)
p.add_argument("--present-fraction", type=float, required=True); p.add_argument("--out", required=True)
a = p.parse_args()
nuc = pd.read_csv(a.nuccore, usecols=["NUCCORE_ACC", "NUCCORE_Description", "NUCCORE_Length", "TAXONOMY_UID"])
tax = pd.read_csv(a.taxonomy, usecols=["TAXONOMY_UID", "TAXONOMY_species"])
nuc = nuc.merge(tax, on="TAXONOMY_UID", how="left").set_index("NUCCORE_ACC")
rows = []
for f in a.screens:
    g = os.path.basename(f).removesuffix(".screen.tsv")
    if os.path.getsize(f) == 0: continue
    s = pd.read_csv(f, sep="\t", header=None, names=["identity", "shared_hashes", "median_multiplicity", "p_value", "plasmid", "comment"])
    for r in s.itertuples():
        n, d = (int(x) for x in r.shared_hashes.split("/"))
        info = nuc.loc[r.plasmid] if r.plasmid in nuc.index else None
        rows.append({"accession": g, "plasmid": r.plasmid, "identity": r.identity, "shared_hashes": r.shared_hashes,
                     "hash_fraction": round(n / d, 3), "p_value": r.p_value, "call": "present" if n / d >= a.present_fraction else "partial",
                     "plasmid_length": None if info is None else int(info.NUCCORE_Length),
                     "plasmid_host": None if info is None else info.TAXONOMY_species,
                     "description": None if info is None else info.NUCCORE_Description})
pd.DataFrame(rows, columns=["accession", "plasmid", "identity", "shared_hashes", "hash_fraction", "p_value", "call",
                            "plasmid_length", "plasmid_host", "description"]).to_csv(a.out, sep="\t", index=False)
print(f"{len(rows)} PLSDB hits")
