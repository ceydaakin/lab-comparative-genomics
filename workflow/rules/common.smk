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
    t = ["results/paper"]
    if config["phylogeny"]["gubbins"]:
        t.append("results/phylogeny/core_fasttree_gubbins.nwk")
    if config["ai_report"]["enabled"]:
        t.append("results/report/results_draft.md")
    return t
