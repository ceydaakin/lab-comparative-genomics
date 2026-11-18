rule download_genomes:
    """Genomes, PGAP proteins and the assembly report from NCBI Datasets."""
    output:
        report="resources/ncbi/assembly_data_report.jsonl",
        fna=expand("resources/genomes/{acc}.fna", acc=ALL),
    params:
        acc=ALL,
    log:
        "logs/download_genomes.log",
    conda:
        "../envs/datasets.yaml"
    shell:
        """
        set -euo pipefail
        tmp=resources/ncbi/tmp && rm -rf $tmp && mkdir -p $tmp resources/ncbi/proteins resources/genomes
        printf '%s\\n' {params.acc} > $tmp/acc.txt
        datasets download genome accession --inputfile $tmp/acc.txt --include genome,protein --filename $tmp/ncbi.zip > {log} 2>&1
        unzip -q -o $tmp/ncbi.zip -d $tmp
        cp $tmp/ncbi_dataset/data/assembly_data_report.jsonl {output.report}
        for a in {params.acc}; do
            cp $tmp/ncbi_dataset/data/$a/*_genomic.fna resources/genomes/$a.fna
            [ -f $tmp/ncbi_dataset/data/$a/protein.faa ] && cp $tmp/ncbi_dataset/data/$a/protein.faa resources/ncbi/proteins/$a.faa || true
        done
        rm -rf $tmp
        """


rule metadata:
    input:
        "resources/ncbi/assembly_data_report.jsonl",
    output:
        "results/metadata.tsv",
    params:
        niches=json.dumps(config["niches"]),
        default=config["default_niche"],
        scripts=SCRIPTS,
    log:
        "logs/metadata.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.scripts}/metadata.py --report {input} --niches '{params.niches}' --default-niche '{params.default}' --out {output} > {log} 2>&1"
