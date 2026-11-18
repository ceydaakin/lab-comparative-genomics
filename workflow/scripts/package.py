"""Paper package: one folder per tool/step (tool name in the folder name), each with its tables and figures/ (PNG + PDF).
Every folder gets a figure; a folder whose figure is missing is reported and the step fails."""
import argparse, glob, os, shutil

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results"); ap.add_argument("--out", required=True)
A = ap.parse_args(); R = A.results

# folder -> (tables/patterns to copy, figure stems)
LAYOUT = {
    "01_QC_BUSCO": (["qc/busco_table.tsv"], ["figS1_busco"]),
    "02_ANI_skani_TYGS": (["ani/ani_matrix.tsv", "ani/clonal_groups.tsv", "tygs/*.tsv", "tygs/*.pdf"], ["fig12_ani_heatmap", "fig16_tygs_dddh"]),
    "03_annotation_Prokka": (["prokka/*/*.txt"], ["fig11_prokka_annotation"]),
    "04_pangenome_Roary_micropan_prabclus": (["roary/summary_statistics.txt", "roary/gene_presence_absence.csv", "roary/gene_presence_absence.Rtab",
                                              "pangenome/micropan/*.tsv"], ["fig2_pan_core", "fig3_pcoa"]),
    "05_phylogeny_FastTree_Gubbins": (["phylogeny/*.nwk", "phylogeny/gubbins.recombination_predictions.gff"], ["fig1_tree_ani"]),
    "06_COG": (["cog/cog_by_set.tsv", "cog/gene_cog.tsv"], ["fig4_cog"]),
    "07_CAZymes_dbCAN": (["cazy/*.tsv"], ["fig6_cazymes"]),
    "08_CRISPR_Cas": (["crispr/crispr_cas.tsv", "crispr/cas_genes.tsv"], ["fig5_crispr_cas"]),
    "09_prophages_PHASTEST": (["phastest/regions.tsv"], ["fig8_prophages"]),
    "10_bacteriocins_BAGEL4": (["bagel4/bagel4_summary.tsv", "bagel4/coverage.txt"], ["fig7_bagel4"]),
    "11_IS_elements_ISEScan": (["isescan/*.tsv"], ["fig9_IS_elements"]),
    "12_AMR_CARD-RGI": (["rgi/rgi_all.tsv"], ["fig14_card_rgi"]),
    "13_plasmids_PLSDB": (["plsdb/plsdb_hits.tsv"], ["fig15_plsdb"]),
    "14_BRIG_comparison": (["brig/focal_unique_regions.tsv"], ["fig10_brig"]),
}
missing, skipped = [], []
for folder, (tables, figs) in LAYOUT.items():
    if not any(glob.glob(f"{R}/{pat}") for pat in tables):
        skipped.append(folder); continue      # step switched off in config (e.g. PHASTEST): no folder
    d = f"{A.out}/{folder}"; os.makedirs(f"{d}/figures", exist_ok=True)
    for pat in tables:
        for f in glob.glob(f"{R}/{pat}"):
            shutil.copy(f, f"{d}/{os.path.basename(f)}")   # file names already carry the accession
    for stem in figs:
        for ext in ("png", "pdf"):
            if os.path.exists(f"{R}/figures/{stem}.{ext}"): shutil.copy(f"{R}/figures/{stem}.{ext}", f"{d}/figures/")
    if not glob.glob(f"{d}/figures/*.png"):
        missing.append(folder)
for f in glob.glob(f"{R}/report/*"): shutil.copy(f, A.out)
if missing:
    raise SystemExit(f"tool folders without a PNG figure: {missing}")
print(f"{len(LAYOUT) - len(skipped)} tool folders in {A.out}; skipped (no outputs): {skipped}")
