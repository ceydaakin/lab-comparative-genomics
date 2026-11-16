rule skani:
    """All-vs-all ANI (stand-in for GET_HOMOLOGUES ANI); species boundary 95-96%."""
    input:
        per_genome("resources/genomes/{acc}.fna"),
    output:
        "results/ani/skani_edges.tsv",
    threads: 8
    log:
        "logs/skani.log",
    conda:
        "../envs/skani.yaml"
    shell:
        "skani triangle -E -t {threads} -o {output} {input} > {log} 2>&1"


rule clonal_groups:
    input:
        edges="results/ani/skani_edges.tsv",
        set="results/qc/genome_set.txt",
    output:
        matrix="results/ani/ani_matrix.tsv",
        groups="results/ani/clonal_groups.tsv",
    params:
        clonal=config["ani"]["clonal_threshold"],
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/clonal_groups.py --skani {input.edges} --genomes {input.set} --clonal {params.clonal} --matrix {output.matrix} --groups {output.groups}"
