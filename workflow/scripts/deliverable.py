"""Hand-off package for collaborators: <Genus>_<species>_results/{00_README.txt, Figures/ (PNG only, descriptive names),
Tables/ (Table 1 + supplementary workbook with a Contents sheet)} and a zip of it. Every number in the README comes from
key_numbers.json, so the text cannot drift from the results."""
import argparse, json, os, shutil
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results"); ap.add_argument("--species", required=True)
ap.add_argument("--zip", required=True, help="path of the zip to write; the unzipped folder sits next to it")
A = ap.parse_args(); R = A.results
k = json.load(open(f"{R}/report/key_numbers.json"))
focal = k["focal"]["strain"]

FIGS = [  # (stem, file name, topic, one-line description)
    ("fig1_tree_ani", "Fig01_Phylogeny_and_ANI", "Species and phylogeny", "Core-genome tree (Roary + FastTree) with ANI heatmap"),
    ("fig12_ani_heatmap", "Fig02_ANI_heatmap", "Species and phylogeny", "Pairwise ANI; box = largest clonal group"),
    ("fig16_tygs_dddh", "Fig03_TYGS_dDDH", "Species and phylogeny", "TYGS digital DNA-DNA hybridisation"),
    ("figS1_busco", "Fig04_BUSCO_quality_control", "Species and phylogeny", "Genome completeness (BUSCO) and the QC cut-off"),
    ("fig2_pan_core", "Fig05_Pangenome_curve", "Pangenome", "Pan/core genome accumulation and Heaps' law (micropan)"),
    ("fig3_pcoa", "Fig06_Gene_content_PCoA", "Pangenome", "Strains by gene content (prabclus Jaccard PCoA)"),
    ("fig4_cog", "Fig07_COG_categories", "Pangenome", "COG categories by pangenome compartment"),
    ("fig11_prokka_annotation", "Fig08_Prokka_annotation", "Pangenome", "Genes, tRNAs and genome size per genome (Prokka)"),
    ("fig10_brig", "Fig09_BRIG_comparison", f"{focal}-specific traits", f"BLASTn rings against {focal} (BRIG method)"),
    ("fig7_bagel4", "Fig10_Bacteriocins_BAGEL4", "Mobile elements and defence", "Bacteriocin areas of interest (BAGEL4)"),
    ("fig5_crispr_cas", "Fig11_CRISPR_Cas", "Mobile elements and defence", "CRISPR spacers and Cas system type"),
    ("fig8_prophages", "Fig12_Prophages", "Mobile elements and defence", "Prophage regions (PHASTEST)"),
    ("fig9_IS_elements", "Fig13_IS_elements", "Mobile elements and defence", "Insertion sequence families (ISEScan)"),
    ("fig15_plsdb", "Fig14_Plasmids_PLSDB", "Mobile elements and defence", "Plasmids detected (PLSDB)"),
    ("fig6_cazymes", "Fig15_CAZymes", "Metabolism and safety", "Carbohydrate-active enzyme families (dbCAN)"),
    ("fig14_card_rgi", "Fig16_Antibiotic_resistance_CARD", "Metabolism and safety", "Antibiotic resistance screen (CARD/RGI)"),
]
SHEETS = {  # sheet label (after 'S<n> ') -> what it shows
    "Genomes & QC": "All input genomes: source, assembly, BUSCO; whether each entered the analysis and why",
    "ANI": "Pairwise ANI (%) matrix", "Roary summary": "Pangenome size by compartment (core, soft-core, shell, cloud)",
    "Focal unique genes": f"Gene families found only in {focal}", "Focal lost genes": f"Families in >= 95% of other strains but absent from {focal}",
    "COG": "COG category share per pangenome compartment (%)", "CRISPR-Cas": "CRISPR arrays, spacers and Cas type per genome",
    "CAZy classes": "CAZyme classes per genome", "CAZy families": "CAZyme families per genome", "BAGEL4": "BAGEL4 bacteriocin areas of interest",
    "Prophages": "PHASTEST prophage regions", "RGI CARD": "CARD/RGI antibiotic-resistance hits", "PLSDB": "PLSDB plasmid matches",
    "IS elements": "ISEScan insertion sequences", "TYGS dDDH": "TYGS dDDH species identification",
}

base = A.zip.removesuffix(".zip")
shutil.rmtree(base, ignore_errors=True); os.makedirs(f"{base}/Figures"); os.makedirs(f"{base}/Tables")
made = []   # numbered over the figures that exist, so a step switched off in config leaves no gap
for stem, name, topic, desc in FIGS:
    src = f"{R}/figures/{stem}.png"
    if os.path.exists(src):
        name = f"Fig{len(made) + 1:02d}_{name.split('_', 1)[1]}"
        shutil.copy(src, f"{base}/Figures/{name}.png"); made.append((name, topic, desc))


def tidy(path, contents=None):
    wb = load_workbook(path)
    if contents is not None:
        ws = wb.create_sheet("Contents", 0); ws.append(["Sheet", "What it shows"])
        for name in wb.sheetnames[1:]:
            ws.append([name, contents.get(name.split(" ", 1)[1], "")])
    for ws in wb.worksheets:
        for c in ws[1]:
            c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="DCE6F2"); c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "A2"
        for i, col in enumerate(ws.iter_cols(min_row=1, max_row=min(ws.max_row, 200)), 1):
            width = max((len(str(c.value)) for c in col if c.value is not None), default=8)
            ws.column_dimensions[get_column_letter(i)].width = min(max(10, width + 2), 60)
    wb.save(path)


shutil.copy(f"{R}/report/Table1_genome_features.xlsx", f"{base}/Tables/Table1_Genome_features.xlsx"); tidy(f"{base}/Tables/Table1_Genome_features.xlsx")
shutil.copy(f"{R}/report/Supplementary_Tables.xlsx", f"{base}/Tables/Supplementary_Tables.xlsx"); tidy(f"{base}/Tables/Supplementary_Tables.xlsx", SHEETS)

q, pan, ani = k["qc"], k["pangenome"], k["ani"]
lines = [f"{A.species.upper()}: POPULATION GENOMICS RESULTS", "", "=" * 70, "1. WHAT WAS DONE", "=" * 70,
         f"- {q['input_genomes']} genomes analysed; {q['analysis_genomes']} passed BUSCO quality control.",
         f"- Focal strain: {focal} ({k['focal']['accession']}), from {k['focal']['niche']}.",
         "- Tools follow the Ortakci-group comparative-genomics pipeline (see docs/tools.md in the repository).",
         "", "=" * 70, "2. FOLDERS", "=" * 70, f"Figures/  {len(made)} figures (PNG, 300 dpi)", "Tables/   Table 1 and the supplementary workbook",
         "", "=" * 70, "3. FIGURES (Figures/)", "=" * 70]
for topic in dict.fromkeys(t for _, t, _ in made):
    lines.append(topic)
    lines += [f"  {n.split('_')[0]}  {d}" for n, t, d in made if t == topic]
    lines.append("")
lines += ["=" * 70, "4. TABLES (Tables/)", "=" * 70,
          "Table1_Genome_features.xlsx    Main-text Table 1: genome features of every strain side by side",
          "Supplementary_Tables.xlsx      Supplementary tables; the first sheet ('Contents') explains each one", "",
          "=" * 70, "5. KEY NUMBERS", "=" * 70,
          f"- One species: all pairwise ANI >= {ani['min_ani']:.1f}%" + (f"; TYGS dDDH {k['tygs']['d4_range_to_type_strain'][0]}-{k['tygs']['d4_range_to_type_strain'][1]}% to the closest type strain" if "tygs" in k else "") + ".",
          f"- Pangenome {pan['Total genes']:,} genes, {pan['Core genes']:,} core; Heaps' alpha {pan['heaps_alpha']} ({'open' if pan['open'] else 'closed'}).",
          f"- {focal}: {k['focal_genes']['unique']} strain-specific gene families (rank {k['focal_genes']['unique_rank']} in the population)."]
if "crispr" in k: lines.append(f"- {focal} CRISPR-Cas: {k['crispr']['focal']['type']}, {k['crispr']['focal']['spacers']} spacers (population median {k['crispr']['spacer_median']:.0f}).")
if "amr" in k: lines.append(f"- Antibiotic resistance: {k['amr']['perfect_hits']} Perfect CARD hits across {k['amr']['genomes']} genomes.")
lines += ["", "The full per-tool folders (tables, PNG and PDF figures, method notes) are in results/paper/."]
open(f"{base}/00_README.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")
shutil.make_archive(base, "zip", os.path.dirname(base) or ".", os.path.basename(base))
print(f"{len(made)} figures, 2 tables -> {base}.zip")
