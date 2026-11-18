rule figures:
    input:
        report_inputs,
    output:
        "results/figures/fig1_tree_ani.png",
    params:
        scripts=SCRIPTS,
        focal=FOCAL,
    log:
        "logs/figures.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/figures.py --results results --focal {params.focal} --out results/figures > {log} 2>&1"


rule tables:
    input:
        report_inputs,
        "results/brig/focal_unique_regions.tsv",
    output:
        xlsx="results/report/Supplementary_Tables.xlsx",
        numbers="results/report/key_numbers.json",
        table1="results/report/Table1_genome_features.xlsx",
    params:
        scripts=SCRIPTS,
        focal=FOCAL,
    log:
        "logs/tables.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/tables.py --results results --focal {params.focal} --xlsx {output.xlsx} --json {output.numbers} --table1 {output.table1} > {log} 2>&1"


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
        scripts=SCRIPTS,
    log:
        "logs/ai_report.log",
    conda:
        "../envs/ai.yaml"
    shell:
        "python {params.scripts}/ai_report.py --numbers {input.numbers} --config '{params.cfg}' --model {params.model} --effort {params.effort} --out {output} > {log} 2>&1"


rule package:
    """One folder per tool (tables + figures/ PNG and PDF) under results/paper; fails if a tool folder has no figure."""
    input:
        "results/figures/fig1_tree_ani.png",
        "results/report/Supplementary_Tables.xlsx",
        "results/brig/focal_unique_regions.tsv",
    output:
        directory("results/paper"),
    params:
        scripts=SCRIPTS,
    log:
        "logs/package.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/package.py --results results --out {output} > {log} 2>&1"


rule deliverable:
    """Hand-off zip for collaborators: 00_README.txt, Figures/ (PNG, descriptive names), Tables/ (Table 1 + supplementary workbook)."""
    input:
        "results/paper",
        "results/report/Table1_genome_features.xlsx",
    output:
        f"{DELIVERABLE}.zip",
    params:
        species=SPECIES,
        scripts=SCRIPTS,
    log:
        "logs/deliverable.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/deliverable.py --results results --species '{params.species}' --zip {output} > {log} 2>&1"
