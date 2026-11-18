rule roary:
    """Roary 3.13.0: -e -n -v, BLASTp identity and core threshold from config (Gumustop & Ortakci 2022)."""
    input:
        per_genome("results/prokka/{acc}/{acc}.gff"),
    output:
        csv="results/roary/gene_presence_absence.csv",
        rtab="results/roary/gene_presence_absence.Rtab",
        summary="results/roary/summary_statistics.txt",
        reference="results/roary/pan_genome_reference.fa",
        core="results/roary/core_gene_alignment.aln",
    params:
        i=config["pangenome"]["roary_identity"],
        cd=config["pangenome"]["roary_core"],
    threads: 8
    log:
        "logs/roary.log",
    conda:
        "../envs/roary.yaml"
    shell:
        # Roary refuses an existing output directory, and Snakemake pre-creates results/roary; run in a temp dir, then move.
        """
        rm -rf results/roary_tmp
        roary -e -n -v -i {params.i} -cd {params.cd} -p {threads} -f results/roary_tmp {input} > {log} 2>&1
        mv results/roary_tmp/* results/roary/ && rm -rf results/roary_tmp
        """


rule micropan_prabclus:
    """Heaps' law (micropan) and Jaccard PCoA (prabclus) on the Roary matrix."""
    input:
        "results/roary/gene_presence_absence.Rtab",
    output:
        heaps="results/pangenome/micropan/heaps.tsv",
        rarefaction="results/pangenome/micropan/rarefaction.tsv",
        pcoa="results/pangenome/micropan/pcoa_jaccard.tsv",
    params:
        nperm=config["pangenome"]["heaps_permutations"],
        scripts=SCRIPTS,
    log:
        "logs/micropan.log",
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript {params.scripts}/micropan_prabclus.R {input} results/pangenome/micropan {params.nperm} > {log} 2>&1"
