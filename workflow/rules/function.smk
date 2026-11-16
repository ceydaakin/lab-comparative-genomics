CDD = "https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd"
COG24 = "https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data"
# dbCAN release used by run_dbcan (dbcan/constants/databases_constants.py, S3_DATABASE_RELEASE)
DBCAN = "https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/dbCAN.hmm"


rule cog_db:
    output:
        pssm=directory("resources/cog/Cog"),
        cddid="resources/cog/cddid.tbl.gz",
        cogdef="resources/cog/cog-24.def.tab",
        cogfun="resources/cog/cog-24.fun.tab",
    conda:
        "../envs/fetch.yaml"
    shell:
        """
        mkdir -p {output.pssm}
        curl -sSfL {CDD}/little_endian/Cog_LE.tar.gz | tar -xz -C {output.pssm}
        curl -sSfL -o {output.cddid} {CDD}/cddid.tbl.gz
        curl -sSfL -o {output.cogdef} {COG24}/cog-24.def.tab
        curl -sSfL -o {output.cogfun} {COG24}/cog-24.fun.tab
        """


rule cog_rpsblast:
    """COG assignment by rpsblast vs NCBI CDD COG PSSMs (stand-in for eggNOG-mapper, whose web form needs an e-mail)."""
    input:
        ref="results/roary/pan_genome_reference.fa",
        db="resources/cog/Cog",
    output:
        faa="results/cog/pan_reference.faa",
        hits="results/cog/rpsblast.tsv",
    threads: 8
    conda:
        "../envs/blast.yaml"
    shell:
        """
        python {SCRIPTS}/translate.py {input.ref} {output.faa}
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
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {SCRIPTS}/cog_annotate.py --rpsblast {input.hits} --reference {input.ref} --presence {input.presence} \
            --cddid {input.cddid} --cog-def {input.cogdef} --cog-fun {input.cogfun} --focal {FOCAL} \
            --genes-out {output.genes} --table-out {output.table}
        """


rule dbcan_db:
    output:
        "resources/dbcan/dbCAN.hmm",
    conda:
        "../envs/fetch.yaml"
    shell:
        "curl -sSfL -o {output} {DBCAN}"


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
    threads: 8
    conda:
        "../envs/python.yaml"
    shell:
        "python {SCRIPTS}/dbcan.py --hmm {input.hmm} --proteomes {input.faa} --evalue {params.e} --coverage {params.c} --threads {threads} --outdir results/cazy"
