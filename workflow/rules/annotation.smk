rule prokka:
    """Prokka 1.14.6, --kingdom Bacteria (all five group papers)."""
    input:
        "resources/genomes/{acc}.fna",
    output:
        gff="results/prokka/{acc}/{acc}.gff",
        faa="results/prokka/{acc}/{acc}.faa",
        txt="results/prokka/{acc}/{acc}.txt",
    params:
        genus=config["species"]["genus"],
        species=config["species"]["species"],
        tag=lambda wc: "L" + "".join(ch for ch in wc.acc if ch.isdigit())[:9],
    threads: 4
    log:
        "logs/prokka/{acc}.log",
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        prokka --kingdom Bacteria --genus {params.genus} --species {params.species} --locustag {params.tag} \
            --prefix {wildcards.acc} --outdir results/prokka/{wildcards.acc} --cpus {threads} --force {input} > {log} 2>&1
        """
