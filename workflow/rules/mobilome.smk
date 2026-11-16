import os


def pgap_proteomes(wildcards):
    """NCBI PGAP proteomes that exist for the analysis set (not every assembly is annotated)."""
    return [f for f in (f"resources/ncbi/proteins/{a}.faa" for a in analysis_set()) if os.path.exists(f)]


# ---------- CRISPR-Cas ----------
rule minced:
    """CRISPR arrays and spacers (MinCED, the CRT-based finder inside Prokka)."""
    input:
        "resources/genomes/{acc}.fna",
    output:
        gff="results/crispr/minced/{acc}.gff",
        txt="results/crispr/minced/{acc}.txt",
    conda:
        "../envs/minced.yaml"
    shell:
        "minced -gffFull {input} {output.txt} {output.gff}"


rule cas_hits:
    input:
        pgap=pgap_proteomes,
        prokka=per_genome("results/prokka/{acc}/{acc}.faa"),
    output:
        ref="results/crispr/pgap_cas.faa",
        query=temp("results/crispr/prokka_all.faa"),
        hits="results/crispr/prokka_vs_pgap_cas.m8",
    threads: 8
    conda:
        "../envs/mmseqs.yaml"
    shell:
        """
        python {SCRIPTS}/cas_inputs.py --pgap {input.pgap} --prokka {input.prokka} --ref-out {output.ref} --query-out {output.query}
        mmseqs easy-search {output.query} {output.ref} {output.hits} results/crispr/tmp --min-seq-id 0.5 -c 0.7 --cov-mode 0 \
            --max-accept 1 --threads {threads} --format-output query,target,pident,qcov > /dev/null
        rm -rf results/crispr/tmp
        """


rule crispr_cas:
    input:
        gff=per_genome("results/crispr/minced/{acc}.gff"),
        hits="results/crispr/prokka_vs_pgap_cas.m8",
        ref="results/crispr/pgap_cas.faa",
    output:
        table="results/crispr/crispr_cas.tsv",
        genes="results/crispr/cas_genes.tsv",
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/crispr_cas.py --minced {input.gff} --hits {input.hits} --cas-ref {input.ref} --out {output.table} --genes-out {output.genes}"


# ---------- Prophages (PHASTEST web API) ----------
rule phastest:
    input:
        "resources/genomes/{acc}.fna",
    output:
        summary="results/phastest/raw/{acc}.summary.txt",
        regions="results/phastest/raw/{acc}.regions.tsv",
    params:
        minutes=config["phastest"]["poll_minutes"],
    resources:
        phastest_jobs=1,        # limit concurrent submissions: snakemake --resources phastest_jobs=4
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/phastest.py --fasta {input} --accession {wildcards.acc} --minutes {params.minutes} --summary {output.summary} --regions {output.regions}"


rule phastest_merge:
    input:
        per_genome("results/phastest/raw/{acc}.regions.tsv"),
    output:
        "results/phastest/regions.tsv",
    shell:
        "awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input} > {output}"


# ---------- Bacteriocins (BAGEL4 web server, no API) ----------
rule bagel4_collect:
    input:
        set="results/qc/genome_set.txt",
    output:
        "results/bagel4/bagel4_summary.tsv",
        "results/bagel4/coverage.txt",
    params:
        sessions=" ".join(config["bagel4"]["sessions"]),
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/bagel4_collect.py --sessions {params.sessions} --genomes {input.set} --outdir results/bagel4"


# ---------- Insertion sequences ----------
rule isescan:
    """ISEScan (ISfinder family HMMs). Single-threaded: multi-threaded FragGeneScan fails silently on some systems."""
    input:
        "resources/genomes/{acc}.fna",
    output:
        "results/isescan/{acc}.tsv",
    log:
        "logs/isescan/{acc}.log",
    conda:
        "../envs/isescan.yaml"
    shell:
        """
        isescan.py --seqfile {input} --output results/isescan/raw/{wildcards.acc} --nthread 1 > {log} 2>&1
        if grep -q "No such file or Empty file" {log}; then echo "ISEScan produced an empty proteome for {wildcards.acc}" >&2; exit 1; fi
        f=$(find results/isescan/raw/{wildcards.acc} -name '*.tsv' | head -1)
        if [ -n "$f" ]; then cp "$f" {output}; else printf 'seqID\\tfamily\\tcluster\\n' > {output}; fi
        """


# ---------- Plasmids (PLSDB) ----------
rule plsdb_db:
    output:
        sketch="resources/plsdb/plsdb_sketch.msh",
        nuccore="resources/plsdb/nuccore.csv",
        taxonomy="resources/plsdb/taxonomy.csv",
    params:
        sketch_url=config["plsdb"]["sketch_url"],
        nuccore_url=config["plsdb"]["nuccore_url"],
        taxonomy_url=config["plsdb"]["taxonomy_url"],
    conda:
        "../envs/fetch.yaml"
    shell:
        """
        curl -sSfL -o {output.sketch} {params.sketch_url}
        curl -sSfL -o {output.nuccore} {params.nuccore_url}
        curl -sSfL -o {output.taxonomy} {params.taxonomy_url}
        """


rule plsdb_screen:
    input:
        fna="resources/genomes/{acc}.fna",
        sketch="resources/plsdb/plsdb_sketch.msh",
    output:
        "results/plsdb/screen/{acc}.screen.tsv",
    params:
        i=config["plsdb"]["min_identity"],
        v=config["plsdb"]["max_pvalue"],
    threads: 4
    conda:
        "../envs/mash.yaml"
    shell:
        "mash screen -w -p {threads} -i {params.i} -v {params.v} {input.sketch} {input.fna} > {output}"


rule plsdb_summary:
    input:
        screens=per_genome("results/plsdb/screen/{acc}.screen.tsv"),
        nuccore="resources/plsdb/nuccore.csv",
        taxonomy="resources/plsdb/taxonomy.csv",
    output:
        "results/plsdb/plsdb_hits.tsv",
    params:
        f=config["plsdb"]["present_hash_fraction"],
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/plsdb_screen.py --screens {input.screens} --nuccore {input.nuccore} --taxonomy {input.taxonomy} --present-fraction {params.f} --out {output}"


# ---------- BRIG-method ring ----------
rule brig:
    input:
        set="results/qc/genome_set.txt",
        tree="results/phylogeny/core_fasttree.nwk",
        meta="results/metadata.tsv",
        fna=per_genome("resources/genomes/{acc}.fna"),
        crispr=f"results/crispr/minced/{FOCAL}.gff",
        prophages=lambda wc: ["results/phastest/regions.tsv"] if config["phastest"]["enabled"] else [],
    output:
        png="results/figures/fig10_brig.png",
        islands="results/brig/focal_unique_regions.tsv",
    params:
        prophages=lambda wc, input: f"--prophages {input.prophages}" if input.prophages else "",
    threads: 8
    conda:
        "../envs/blast.yaml"
    shell:
        """
        python {SCRIPTS}/brig.py --focal {FOCAL} --genome-dir resources/genomes --genomes {input.set} --tree {input.tree} \
            --metadata {input.meta} --workdir results/brig/blast {params.prophages} --crispr-gff {input.crispr} \
            --threads {threads} --figure results/figures/fig10_brig --islands {output.islands}
        """
