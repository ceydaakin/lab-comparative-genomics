rule cog_db:
    output:
        pssm=directory("resources/cog/Cog"),
        cddid="resources/cog/cddid.tbl.gz",
        cogdef="resources/cog/cog-24.def.tab",
        cogfun="resources/cog/cog-24.fun.tab",
    params:
        cdd=CDD,
        cog24=COG24,
    log:
        "logs/cog_db.log",
    conda:
        "../envs/fetch.yaml"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output.pssm}
        curl -sSfL {params.cdd}/little_endian/Cog_LE.tar.gz | tar -xz -C {output.pssm}
        curl -sSfL -o {output.cddid} {params.cdd}/cddid.tbl.gz
        curl -sSfL -o {output.cogdef} {params.cog24}/cog-24.def.tab
        curl -sSfL -o {output.cogfun} {params.cog24}/cog-24.fun.tab
        """


rule cog_rpsblast:
    """COG assignment by rpsblast vs NCBI CDD COG PSSMs (stand-in for eggNOG-mapper, whose web form needs an e-mail)."""
    input:
        ref="results/roary/pan_genome_reference.fa",
        db="resources/cog/Cog",
    output:
        faa="results/cog/pan_reference.faa",
        hits="results/cog/rpsblast.tsv",
    params:
        scripts=SCRIPTS,
    log:
        "logs/cog_rpsblast.log",
    threads: 8
    conda:
        "../envs/blast.yaml"
    shell:
        """
        exec > {log} 2>&1
        python {params.scripts}/translate.py {input.ref} {output.faa}
        rpsblast -query {output.faa} -db {input.db}/Cog -evalue 1e-5 -max_target_seqs 1 \
            -outfmt '6 qseqid sseqid evalue bitscore' -num_threads {threads} -out {output.hits}
        """


rule cog_annotate:
    input:
        hits="results/cog/rpsblast.tsv",
        ref="results/roary/pan_genome_reference.fa",
        presence="results/roary/gene_presence_absence.csv",
        cddid="resources/cog/cddid.tbl.gz",
        cogdef="resources/cog/cog-24.def.tab",
        cogfun="resources/cog/cog-24.fun.tab",
    output:
        genes="results/cog/gene_cog.tsv",
        table="results/cog/cog_by_set.tsv",
    params:
        scripts=SCRIPTS,
        focal=FOCAL,
    log:
        "logs/cog_annotate.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        exec > {log} 2>&1
        python {params.scripts}/cog_annotate.py --rpsblast {input.hits} --reference {input.ref} --presence {input.presence} \
            --cddid {input.cddid} --cog-def {input.cogdef} --cog-fun {input.cogfun} --focal {params.focal} \
            --genes-out {output.genes} --table-out {output.table}
        """


rule dbcan_db:
    output:
        "resources/dbcan/dbCAN.hmm",
    params:
        dbcan_url=DBCAN,
    log:
        "logs/dbcan_db.log",
    conda:
        "../envs/fetch.yaml"
    shell:
        "curl -sSfL -o {output} {params.dbcan_url} > {log} 2>&1"


rule dbcan:
    """CAZymes with dbCAN HMMs on Prokka proteins (Gumustop & Ortakci 2022 thresholds)."""
    input:
        hmm="resources/dbcan/dbCAN.hmm",
        faa=per_genome("results/prokka/{acc}/{acc}.faa"),
    output:
        "results/cazy/cazy_family_counts.tsv",
        "results/cazy/cazy_class_counts.tsv",
        "results/cazy/cazy_hits.tsv",
    params:
        e=config["cazy"]["evalue"],
        c=config["cazy"]["coverage"],
        scripts=SCRIPTS,
    log:
        "logs/dbcan.log",
    threads: 8
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/dbcan.py --hmm {input.hmm} --proteomes {input.faa} --evalue {params.e} --coverage {params.c} --threads {threads} --outdir results/cazy > {log} 2>&1"
