"""Supplementary_Tables.xlsx (one sheet per table) and key_numbers.json (the computed facts the AI report may use).
Every number the report can cite comes from this file; the AI step never sees raw data it could misread."""
import argparse, glob, json, os
import pandas as pd
from scipy.stats import mannwhitneyu

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results"); ap.add_argument("--focal", required=True)
ap.add_argument("--xlsx", required=True); ap.add_argument("--json", required=True)
A = ap.parse_args(); R = A.results; F = A.focal
G = [l.strip() for l in open(f"{R}/qc/genome_set.txt") if l.strip()]
meta = pd.read_csv(f"{R}/metadata.tsv", sep="\t").set_index("accession")
strain = lambda a: meta.strain.get(a, a)
opt = lambda p: pd.read_csv(p, sep="\t") if os.path.exists(p) and os.path.getsize(p) else None
sheets, k = {}, {"focal": {"accession": F, "strain": strain(F), "niche": meta.niche.get(F)}}

busco = pd.read_csv(f"{R}/qc/busco_table.tsv", sep="\t"); sheets["S1 Genomes & QC"] = busco
k["qc"] = {"input_genomes": len(busco), "analysis_genomes": len(G), "niches": meta.niche.reindex(G).value_counts().to_dict()}
ani = pd.read_csv(f"{R}/ani/ani_matrix.tsv", sep="\t", index_col=0); sheets["S2 ANI"] = ani.reset_index()
cg = pd.read_csv(f"{R}/ani/clonal_groups.tsv", sep="\t"); vc = cg.clonal_group.value_counts()
big = cg[cg.clonal_group == vc.index[0]].accession.tolist()
clonal = set(big) if len(big) >= 3 else set()
k["ani"] = {"min_ani": float(ani.where(ani < 100).min().min()), "focal_ani_range": [float(ani.loc[F].drop(F).min()), float(ani.loc[F].drop(F).max())],
            "clonal_groups": int(cg.clonal_group.nunique()), "largest_clonal_group": [strain(a) for a in big]}
st = pd.read_csv(f"{R}/roary/summary_statistics.txt", sep="\t", header=None, names=["category", "definition", "genes"]); sheets["S3 Roary summary"] = st
h = pd.read_csv(f"{R}/pangenome/micropan/heaps.tsv", sep="\t").set_index("parameter").value
pc = pd.read_csv(f"{R}/pangenome/micropan/pcoa_jaccard.tsv", sep="\t")
k["pangenome"] = {**{r.category: int(r.genes) for r in st.itertuples()}, "heaps_alpha": round(float(h["alpha"]), 3),
                  "open": bool(h["alpha"] < 1), "pcoa_variance": [round(float(pc.var1.iloc[0]), 1), round(float(pc.var2.iloc[0]), 1)]}
pa = pd.read_csv(f"{R}/roary/gene_presence_absence.csv", low_memory=False)
gcols = pa.columns[pa.columns.get_loc("Avg group size nuc") + 1:]; fc = [c for c in gcols if c.startswith(F)][0]; others = [c for c in gcols if c != fc]
pres = pa[gcols].notna()
uniq = pa[pres[fc] & (pres[others].sum(axis=1) == 0)][["Gene", "Annotation"]]
lost = pa[~pres[fc] & (pres[others].mean(axis=1) >= 0.95)][["Gene", "Annotation", "No. isolates"]]
sheets["S4 Focal unique genes"], sheets["S5 Focal lost genes"] = uniq, lost
per_genome_unique = pres[pres.sum(axis=1) == 1].sum().sort_values(ascending=False)
k["focal_genes"] = {"unique": len(uniq), "unique_rank": int(list(per_genome_unique.index).index(fc)) + 1,
                    "unique_annotations_top": uniq.Annotation.value_counts().head(15).to_dict(), "lost_annotations": lost.Annotation.tolist()}
cog = opt(f"{R}/cog/cog_by_set.tsv")
if cog is not None:
    sheets["S6 COG"] = cog
    k["cog_percent"] = cog[cog.genes > 0].pivot(index="cog_category", columns="set", values="percent_of_set").fillna(0).round(1).to_dict()
cr = opt(f"{R}/crispr/crispr_cas.tsv")
if cr is not None:
    sheets["S7 CRISPR-Cas"] = cr; f = cr.set_index("accession").loc[F]
    k["crispr"] = {"type_counts": cr.cas_type.value_counts().to_dict(), "spacer_median": float(cr.spacers.median()),
                   "focal": {"type": f.cas_type, "arrays": int(f.arrays), "spacers": int(f.spacers)}}
cz = opt(f"{R}/cazy/cazy_class_counts.tsv")
if cz is not None:
    cz = cz.set_index(cz.columns[0]); sheets["S8 CAZy classes"] = cz.reset_index(); tot = cz.sum(axis=1)
    fam = pd.read_csv(f"{R}/cazy/cazy_family_counts.tsv", sep="\t", index_col=0); sheets["S9 CAZy families"] = fam.reset_index()
    only = fam.index[(fam.get(F, 0) > 0) & (fam.drop(columns=[F], errors="ignore").sum(axis=1) == 0)].tolist()
    k["cazy"] = {"total_range": [int(tot.min()), int(tot.max())], "median": float(tot.median()), "focal": cz.loc[F].astype(int).to_dict(), "focal_only_families": only}
bg = opt(f"{R}/bagel4/bagel4_summary.tsv")
if bg is not None:
    sheets["S10 BAGEL4"] = bg
    k["bagel4"] = {"coverage": open(f"{R}/bagel4/coverage.txt").read().strip(), "aois": bg.assign(strain=bg.accession.map(strain))[["strain", "region_size_bp", "best_hit"]].to_dict("records")}
ph = opt(f"{R}/phastest/regions.tsv")
if ph is not None:
    ph = ph[ph.accession.isin(G)]; sheets["S11 Prophages"] = ph
    it = ph[ph.completeness == "intact"].groupby("accession").size().reindex(G, fill_value=0)
    k["prophages"] = {"regions": len(ph), "intact": int(it.sum()), "genomes_with_intact": int((it > 0).sum()),
                      "focal_intact": ph[(ph.accession == F) & (ph.completeness == "intact")][["length_kb", "position"]].to_dict("records")}
    if clonal and len(clonal) < len(G):
        a, b = it[[g for g in G if g in clonal]], it[[g for g in G if g not in clonal]]
        u = mannwhitneyu(a, b)
        k["prophages"]["clonal_vs_rest"] = {"clonal_with_intact": f"{(a > 0).sum()}/{len(a)}", "rest_with_intact": f"{(b > 0).sum()}/{len(b)}",
                                            "mann_whitney_U": float(u.statistic), "p": float(u.pvalue)}
rgi = opt(f"{R}/rgi/rgi_all.tsv")
if rgi is not None:
    sheets["S12 RGI CARD"] = rgi
    k["amr"] = {"perfect_hits": int((rgi.cut_off == "Perfect").sum()), "strict_hits": rgi[rgi.cut_off == "Strict"].groupby("best_hit").accession.nunique().to_dict(),
                "genomes": len(G)}
pl = opt(f"{R}/plsdb/plsdb_hits.tsv")
if pl is not None:
    sheets["S13 PLSDB"] = pl
    k["plasmids"] = {"focal_present": pl[(pl.accession == F) & (pl.call == "present")][["plasmid", "hash_fraction", "plasmid_length", "plasmid_host"]].to_dict("records"),
                     "hosts_present": pl[pl.call == "present"].plasmid_host.value_counts().to_dict()}
iss = [pd.read_csv(f, sep="\t").assign(accession=os.path.basename(f)[:-4]) for f in sorted(glob.glob(f"{R}/isescan/*.tsv")) if os.path.getsize(f)]
if iss:
    isd = pd.concat(iss); sheets["S14 IS elements"] = isd; n = isd.groupby("accession").size().reindex(G, fill_value=0)
    k["is_elements"] = {"median_per_genome": float(n.median()), "range": [int(n.min()), int(n.max())], "focal_total": int(n.get(F, 0)),
                        "focal_families": isd[isd.accession == F].family.value_counts().to_dict()}
isl = opt(f"{R}/brig/focal_unique_regions.tsv")
if isl is not None: k["brig_focal_only_regions"] = isl.to_dict("records")

with pd.ExcelWriter(A.xlsx) as xw:
    for name, df in sheets.items(): df.to_excel(xw, sheet_name=name[:31], index=False)
json.dump(k, open(A.json, "w"), indent=1, default=str)
print(f"{len(sheets)} supplementary sheets; key numbers -> {A.json}")
