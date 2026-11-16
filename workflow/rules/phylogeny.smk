rule fasttree_core:
    """Roary core-gene alignment -> FastTree 2 GTR (Gumustop & Ortakci 2022, 2023)."""
    input:
        "results/roary/core_gene_alignment.aln",
    output:
        "results/phylogeny/core_fasttree.nwk",
    log:
        "logs/fasttree.log",
    conda:
        "../envs/phylo.yaml"
    shell:
        "FastTree -gtr -nt {input} > {output} 2> {log}"


rule gubbins:
    """Recombination filtering (Gubbins), then FastTree on the filtered polymorphic sites (Gumustop et al. 2025; Kurt et al. 2026).
    The bioconda Gubbins build segfaults on Apple Silicon: run this on Linux."""
    input:
        "results/roary/core_gene_alignment.aln",
    output:
        tree="results/phylogeny/core_fasttree_gubbins.nwk",
        gff="results/phylogeny/gubbins.recombination_predictions.gff",
    threads: 4
    log:
        "logs/gubbins.log",
    conda:
        "../envs/phylo.yaml"
    shell:
        """
        cd results/phylogeny
        run_gubbins.py --prefix gubbins --tree-builder fasttree --threads {threads} ../../{input} > ../../{log} 2>&1
        FastTree -gtr -nt gubbins.filtered_polymorphic_sites.fasta > core_fasttree_gubbins.nwk 2>> ../../{log}
        """
