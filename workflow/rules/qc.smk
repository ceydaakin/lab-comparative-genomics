rule busco_lineage:
    output:
        directory(f"resources/busco/lineages/{LINEAGE}"),
    log:
        "logs/busco_lineage.log",
    conda:
        "../envs/busco.yaml"
    shell:
        f"busco --download {LINEAGE} --download_path resources/busco > {{log}} 2>&1"


rule busco:
    """BUSCO completeness, genome mode (Gumustop & Ortakci 2022; Kurt et al. 2026)."""
    input:
        fna="resources/genomes/{acc}.fna",
        lineage=rules.busco_lineage.output,
    output:
        "results/qc/busco/{acc}.json",
    threads: 2
    log:
        "logs/busco/{acc}.log",
    conda:
        "../envs/busco.yaml"
    shell:
        """
        busco -i {input.fna} -m genome -l {input.lineage} --offline -o {wildcards.acc} --out_path results/qc/busco/run -c {threads} -f > {log} 2>&1
        cp results/qc/busco/run/{wildcards.acc}/short_summary.specific.*.{wildcards.acc}.json {output}
        """


checkpoint qc_select:
    input:
        jsons=expand("results/qc/busco/{acc}.json", acc=ALL),
        meta="results/metadata.tsv",
    output:
        table="results/qc/busco_table.tsv",
        set="results/qc/genome_set.txt",
        cut="results/qc/min_busco.txt",
    params:
        cut=config["qc"]["min_busco_complete"],
        exclude=" ".join(config["qc"].get("exclude", [])),
        scripts=SCRIPTS,
        focal=FOCAL,
    log:
        "logs/qc_select.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        exec > {log} 2>&1
        python {params.scripts}/busco_select.py --summaries {input.jsons} --metadata {input.meta} --min-complete {params.cut} \
            --exclude {params.exclude} --focal {params.focal} --table {output.table} --set {output.set}
        echo {params.cut} > {output.cut}
        """
