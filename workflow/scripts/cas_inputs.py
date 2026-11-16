"""Inputs for Cas typing: (1) Cas-named proteins from the NCBI PGAP proteomes of the input genomes (the reference set),
(2) all Prokka proteins of the analysis set, ids prefixed '<accession>|' so hits map back to genomes."""
import argparse, os

p = argparse.ArgumentParser()
p.add_argument("--pgap", nargs="*", default=[]); p.add_argument("--prokka", nargs="+", required=True)
p.add_argument("--ref-out", required=True); p.add_argument("--query-out", required=True)
a = p.parse_args()


def is_cas(name):
    low = name.lower().removeprefix("mag: ")
    return "crispr" in low or low.startswith("cas") or "csn2" in low or "cst" in low or "csm" in low


n = 0
with open(a.ref_out, "w") as o:
    for f in a.pgap:
        keep = False
        for line in open(f):
            if line.startswith(">"):
                pid, _, rest = line[1:].strip().partition(" ")
                name = rest.split(" [")[0]
                keep = is_cas(name)
                if keep: o.write(f">{pid} {name}\n"); n += 1
            elif keep:
                o.write(line)
with open(a.query_out, "w") as o:
    for f in a.prokka:
        g = os.path.basename(f).removesuffix(".faa")
        for line in open(f):
            o.write(f">{g}|{line[1:].split()[0]}\n" if line.startswith(">") else line)
if n == 0:
    raise SystemExit("no Cas-named PGAP proteins found: Cas typing needs at least one NCBI-annotated genome with Cas genes")
print(f"{n} Cas reference proteins")
