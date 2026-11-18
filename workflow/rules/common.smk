import glob
import json
import os
import pandas as pd

ALL = pd.read_csv(config["genomes"], sep="\t").accession.dropna().astype(str).tolist()
FOCAL = config["focal"]
if FOCAL not in ALL:
    raise ValueError(f"focal genome {FOCAL} is not listed in {config['genomes']}")
SCRIPTS = os.path.join(workflow.basedir, "scripts")  # workflow.basedir = folder of the Snakefile


def analysis_set(wildcards=None):
    """Genomes passing BUSCO QC; resolved after the qc_select checkpoint has run."""
    with open(checkpoints.qc_select.get().output.set) as f:
        return [l.strip() for l in f if l.strip()]


def per_genome(pattern):
    return lambda wc: expand(pattern, acc=analysis_set())


def final_targets(wildcards):
    t = ["results/paper", f"{DELIVERABLE}.zip"]
    if config["phylogeny"]["gubbins"]:
        t.append("results/phylogeny/core_fasttree_gubbins.nwk")
    if config["ai_report"]["enabled"]:
        t.append("results/report/results_draft.md")
    return t


# ---- shared constants and input functions; the other .smk files contain only workflow steps (snakemake --lint) ----
CDD = "https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd"
COG24 = "https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data"
# dbCAN release used by run_dbcan (dbcan/constants/databases_constants.py, S3_DATABASE_RELEASE)
DBCAN = "https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/dbCAN.hmm"


LINEAGE = config["qc"]["busco_lineage"]


def pgap_proteomes(wildcards):
    """NCBI PGAP proteomes that exist for the analysis set (not every assembly is annotated)."""
    return [f for f in (f"resources/ncbi/proteins/{a}.faa" for a in analysis_set()) if os.path.exists(f)]


def report_inputs(wildcards):
    i = ["results/qc/busco_table.tsv", "results/ani/ani_matrix.tsv", "results/ani/clonal_groups.tsv", "results/phylogeny/core_fasttree.nwk",
         "results/pangenome/micropan/heaps.tsv", "results/roary/gene_presence_absence.csv", "results/cog/cog_by_set.tsv",
         "results/crispr/crispr_cas.tsv", "results/cazy/cazy_family_counts.tsv", "results/rgi/rgi_all.tsv",
         "results/plsdb/plsdb_hits.tsv", "results/bagel4/bagel4_summary.tsv", "results/tygs/tygs_species_identification.tsv"]
    i += expand("results/isescan/{acc}.tsv", acc=analysis_set())
    if config["phastest"]["enabled"]:
        i.append("results/phastest/regions.tsv")
    return i


TOOLS = ["BUSCO 6.1", "skani", "Prokka 1.14.6", "Roary 3.13.0", "micropan", "prabclus", "FastTree 2", "Gubbins", "rpsblast vs NCBI COG",
         "dbCAN HMMs", "MinCED", "PHASTEST", "BAGEL4", "ISEScan", "RGI/CARD", "PLSDB (mash screen)", "BLASTn BRIG-method rings"]


SPECIES = f"{config['species']['genus']} {config['species']['species']}"
DELIVERABLE = f"results/deliverable/{config['species']['genus']}_{config['species']['species']}_results"

