"""CAZyme annotation with the dbCAN HMM database (same HMMs as the dbCAN3 server) via pyhmmer.
Thresholds as in Gumustop & Ortakci 2022: domain E-value and HMM coverage cut-offs; best family per protein."""
import argparse, itertools, os
import pandas as pd, pyhmmer

p = argparse.ArgumentParser()
p.add_argument("--hmm", required=True); p.add_argument("--proteomes", nargs="+", required=True)
p.add_argument("--evalue", type=float, required=True); p.add_argument("--coverage", type=float, required=True)
p.add_argument("--threads", type=int, default=4); p.add_argument("--outdir", required=True)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)
with pyhmmer.plan7.HMMFile(a.hmm) as f: hmms = list(f)
alpha = pyhmmer.easel.Alphabet.amino()
rows = []
for faa in a.proteomes:
    g = os.path.basename(faa).removesuffix(".faa")
    with pyhmmer.easel.SequenceFile(faa, digital=True, alphabet=alpha) as sf: seqs = sf.read_block()
    best = {}
    for hits in pyhmmer.hmmsearch(hmms, seqs, cpus=a.threads, E=a.evalue):
        fam = hits.query.name if isinstance(hits.query.name, str) else hits.query.name.decode()
        fam = fam.replace(".hmm", "")
        for h in hits:
            name = h.name if isinstance(h.name, str) else h.name.decode()
            for d in h.domains:
                cov = (d.alignment.hmm_to - d.alignment.hmm_from + 1) / hits.query.M
                if d.i_evalue < a.evalue and cov > a.coverage and (name not in best or d.i_evalue < best[name][1]):
                    best[name] = (fam, d.i_evalue)
    rows += [(g, prot, fam) for prot, (fam, _) in best.items()]
df = pd.DataFrame(rows, columns=["genome", "protein", "family"])
df.to_csv(f"{a.outdir}/cazy_hits.tsv", sep="\t", index=False)
df["cls"] = df.family.map(lambda f: "".join(itertools.takewhile(str.isalpha, f)))  # class prefix: GH, GT, CE, PL, AA, CBM
pd.crosstab(df.family, df.genome).to_csv(f"{a.outdir}/cazy_family_counts.tsv", sep="\t")
pd.crosstab(df.genome, df.cls).to_csv(f"{a.outdir}/cazy_class_counts.tsv", sep="\t")
print(f"{len(df)} CAZyme domains in {df.genome.nunique()} genomes")
