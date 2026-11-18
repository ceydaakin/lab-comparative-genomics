rule card_db:
    """CARD (latest) loaded as an RGI local database in resources/card/localDB."""
    output:
        directory("resources/card/localDB"),
    log:
        "logs/card_db.log",
    conda:
        "../envs/rgi.yaml"
    shell:
        """
        cd resources/card
        curl -sSfL -o card.tar.bz2 https://card.mcmaster.ca/latest/data && tar -xjf card.tar.bz2 ./card.json && rm card.tar.bz2
        rgi load --card_json card.json --local > ../../{log} 2>&1
        """


rule rgi:
    """RGI against CARD, Strict/Perfect hits (Gumustop & Ortakci 2022, 2023)."""
    input:
        faa="results/prokka/{acc}/{acc}.faa",
        db="resources/card/localDB",
    output:
        "results/rgi/{acc}.txt",
    threads: 4
    log:
        "logs/rgi/{acc}.log",
    conda:
        "../envs/rgi.yaml"
    shell:
        # --local reads ./localDB, so run from resources/card
        """
        cd resources/card
        rgi main -i ../../{input.faa} -o ../../results/rgi/{wildcards.acc} -t protein -a DIAMOND -n {threads} --local --clean > ../../{log} 2>&1
        """


rule rgi_merge:
    input:
        per_genome("results/rgi/{acc}.txt"),
    output:
        "results/rgi/rgi_all.tsv",
    params:
        scripts=SCRIPTS,
    log:
        "logs/rgi_merge.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/rgi_merge.py --reports {input} --out {output} > {log} 2>&1"
