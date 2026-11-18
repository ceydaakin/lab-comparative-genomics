# ---------- CRISPR-Cas ----------
rule minced:
    """CRISPR arrays and spacers (MinCED, the CRT-based finder inside Prokka)."""
    input:
        "resources/genomes/{acc}.fna",
    output:
        gff="results/crispr/minced/{acc}.gff",
        txt="results/crispr/minced/{acc}.txt",
    log:
        "logs/minced/{acc}.log",
    conda:
        "../envs/minced.yaml"
    shell:
        "minced -gffFull {input} {output.txt} {output.gff} > {log} 2>&1"


rule cas_hits:
    input:
        pgap=pgap_proteomes,
        prokka=per_genome("results/prokka/{acc}/{acc}.faa"),
    output:
        ref="results/crispr/pgap_cas.faa",
        query=temp("results/crispr/prokka_all.faa"),
        hits="results/crispr/prokka_vs_pgap_cas.m8",
    params:
        scripts=SCRIPTS,
    log:
        "logs/cas_hits.log",
    threads: 8
    conda:
        "../envs/mmseqs.yaml"
    shell:
        """
        exec > {log} 2>&1
        python {params.scripts}/cas_inputs.py --pgap {input.pgap} --prokka {input.prokka} --ref-out {output.ref} --query-out {output.query}
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
    params:
        scripts=SCRIPTS,
    log:
        "logs/crispr_cas.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/crispr_cas.py --minced {input.gff} --hits {input.hits} --cas-ref {input.ref} --out {output.table} --genes-out {output.genes} > {log} 2>&1"


# ---------- Prophages (PHASTEST web API) ----------
rule phastest:
    input:
        "resources/genomes/{acc}.fna",
    output:
        summary="results/phastest/raw/{acc}.summary.txt",
        regions="results/phastest/raw/{acc}.regions.tsv",
    params:
        minutes=config["phastest"]["poll_minutes"],
        scripts=SCRIPTS,
    log:
        "logs/phastest/{acc}.log",
    resources:
        phastest_jobs=1,        # limit concurrent submissions: snakemake --resources phastest_jobs=4
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/phastest.py --fasta {input} --accession {wildcards.acc} --minutes {params.minutes} --summary {output.summary} --regions {output.regions} > {log} 2>&1"


rule phastest_merge:
    input:
        per_genome("results/phastest/raw/{acc}.regions.tsv"),
    output:
        "results/phastest/regions.tsv",
    log:
        "logs/phastest_merge.log",
    conda:
        "../envs/python.yaml"
    shell:
        "awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input} > {output} 2> {log}"


# ---------- Bacteriocins (BAGEL4 web server, no API) ----------
rule bagel4_collect:
    input:
        set="results/qc/genome_set.txt",
    output:
        "results/bagel4/bagel4_summary.tsv",
        "results/bagel4/coverage.txt",
    params:
        sessions=" ".join(config["bagel4"]["sessions"]),
        scripts=SCRIPTS,
    log:
        "logs/bagel4_collect.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/bagel4_collect.py --sessions {params.sessions} --genomes {input.set} --outdir results/bagel4 > {log} 2>&1"


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
    log:
        "logs/plsdb_db.log",
    conda:
        "../envs/fetch.yaml"
    shell:
        """
        exec > {log} 2>&1
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
    log:
        "logs/plsdb_screen/{acc}.log",
    threads: 4
    conda:
        "../envs/mash.yaml"
    shell:
        "mash screen -w -p {threads} -i {params.i} -v {params.v} {input.sketch} {input.fna} > {output} 2> {log}"


rule plsdb_summary:
    input:
        screens=per_genome("results/plsdb/screen/{acc}.screen.tsv"),
        nuccore="resources/plsdb/nuccore.csv",
        taxonomy="resources/plsdb/taxonomy.csv",
    output:
        "results/plsdb/plsdb_hits.tsv",
    params:
        f=config["plsdb"]["present_hash_fraction"],
        scripts=SCRIPTS,
    log:
        "logs/plsdb_summary.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/plsdb_screen.py --screens {input.screens} --nuccore {input.nuccore} --taxonomy {input.taxonomy} --present-fraction {params.f} --out {output} > {log} 2>&1"


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
        scripts=SCRIPTS,
        focal=FOCAL,
    log:
        "logs/brig.log",
    threads: 8
    conda:
        "../envs/blast.yaml"
    shell:
        """
        exec > {log} 2>&1
        python {params.scripts}/brig.py --focal {params.focal} --genome-dir resources/genomes --genomes {input.set} --tree {input.tree} \
            --metadata {input.meta} --workdir results/brig/blast {params.prophages} --crispr-gff {input.crispr} \
            --threads {threads} --figure results/figures/fig10_brig --islands {output.islands}
        """


# ---------- Species: TYGS dDDH (web form, collected by job GUID) ----------
rule tygs_collect:
    input:
        set="results/qc/genome_set.txt",
        meta="results/metadata.tsv",
    output:
        "results/tygs/tygs_species_identification.tsv",
        "results/tygs/coverage.txt",
    params:
        jobs=" ".join(config.get("tygs", {}).get("jobs", [])),
        scripts=SCRIPTS,
    log:
        "logs/tygs_collect.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/tygs_collect.py --jobs {params.jobs} --genomes {input.set} --metadata {input.meta} --outdir results/tygs > {log} 2>&1"
