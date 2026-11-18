import glob


def report_inputs(wildcards):
    i = ["results/qc/busco_table.tsv", "results/ani/ani_matrix.tsv", "results/ani/clonal_groups.tsv", "results/phylogeny/core_fasttree.nwk",
         "results/pangenome/micropan/heaps.tsv", "results/roary/gene_presence_absence.csv", "results/cog/cog_by_set.tsv",
         "results/crispr/crispr_cas.tsv", "results/cazy/cazy_family_counts.tsv", "results/rgi/rgi_all.tsv",
         "results/plsdb/plsdb_hits.tsv", "results/bagel4/bagel4_summary.tsv"]
    i += expand("results/isescan/{acc}.tsv", acc=analysis_set())
    if config["phastest"]["enabled"]:
        i.append("results/phastest/regions.tsv")
    return i


rule figures:
    input:
        report_inputs,
    output:
        "results/figures/fig1_tree_ani.png",
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/figures.py --results results --focal {FOCAL} --out results/figures"


rule tables:
    input:
        report_inputs,
        "results/brig/focal_unique_regions.tsv",
    output:
        xlsx="results/report/Supplementary_Tables.xlsx",
        numbers="results/report/key_numbers.json",
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/tables.py --results results --focal {FOCAL} --xlsx {output.xlsx} --json {output.numbers}"


TOOLS = ["BUSCO 6.1", "skani", "Prokka 1.14.6", "Roary 3.13.0", "micropan", "prabclus", "FastTree 2", "Gubbins", "rpsblast vs NCBI COG",
         "dbCAN HMMs", "MinCED", "PHASTEST", "BAGEL4", "ISEScan", "RGI/CARD", "PLSDB (mash screen)", "BLASTn BRIG-method rings"]


rule ai_report:
    """Optional: Claude drafts Results + figure legends from key_numbers.json (needs Anthropic credentials)."""
    input:
        numbers="results/report/key_numbers.json",
        fig="results/figures/fig1_tree_ani.png",
    output:
        "results/report/results_draft.md",
    params:
        cfg=lambda wc: json.dumps({"species": f"{config['species']['genus']} {config['species']['species']}", "tools": TOOLS,
                                   "figures": sorted(os.path.basename(f)[:-4] for f in glob.glob("results/figures/*.png"))}),
        model=config["ai_report"]["model"],
        effort=config["ai_report"]["effort"],
    conda:
        "../envs/ai.yaml"
    shell:
        "python {SCRIPTS}/ai_report.py --numbers {input.numbers} --config '{params.cfg}' --model {params.model} --effort {params.effort} --out {output}"


rule package:
    """One folder per tool (tables + figures/ PNG and PDF) under results/paper; fails if a tool folder has no figure."""
    input:
        "results/figures/fig1_tree_ani.png",
        "results/report/Supplementary_Tables.xlsx",
        "results/brig/focal_unique_regions.tsv",
    output:
        directory("results/paper"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/package.py --results results --out {output}"
