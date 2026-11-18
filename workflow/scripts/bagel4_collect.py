"""Collect BAGEL4 web-server sessions (bagel4.molgenrug.nl has no API) and summarise areas of interest (AOIs) per genome.
Core peptide = an AOI GeneTable row BAGEL4 colours as core (green) or that carries a scored database hit."""
import argparse, os, urllib.request
from html.parser import HTMLParser
import pandas as pd
from parsers import acc_from_upload, db_hits

ap = argparse.ArgumentParser()
ap.add_argument("--sessions", nargs="*", default=[]); ap.add_argument("--genomes", required=True); ap.add_argument("--outdir", required=True)
args = ap.parse_args()
OUT = args.outdir
os.makedirs(f"{OUT}/raw", exist_ok=True)
G = [l.strip() for l in open(args.genomes) if l.strip()]

class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.files = []
    def handle_starttag(self, tag, attrs):
        href = dict(attrs).get("href", "")
        if tag == "a" and href and not href.startswith(("?", "/")):
            self.files.append(href)




rows, maps = [], []
for session in args.sessions:
    base = f"https://bagel4.molgenrug.nl/results/{session}/"
    get = lambda name: urllib.request.urlopen(base + name, timeout=120).read().decode("utf-8", "replace")
    p = Links(); p.feed(get("")); files = p.files
    if "00.OverviewGeneTables.json" not in files:
        raise SystemExit(f"BAGEL4 session {session} has not finished (no 00.OverviewGeneTables.json yet)")
    contigs = pd.read_csv(base + "00.all_contigs.table", sep="\t")
    contigs["accession"] = contigs.filename.map(acc_from_upload); contigs["session"] = session; maps.append(contigs)
    q2acc = dict(zip(contigs.queryname, contigs.accession))
    for f in [f for f in files if f.endswith(".GeneTable") and "AOI_" in f]:
        open(f"{OUT}/raw/{f}", "w").write(get(f))
        t = pd.read_csv(f"{OUT}/raw/{f}", sep="\t", dtype=str).fillna("")
        query, aoi = f.rsplit(".AOI_", 1)
        core = t[t.gene_color.isin(["green", "#70c947"]) | t.function.map(lambda fn: bool(db_hits(fn)))]
        hits = sorted({h for fn in core.function for h in db_hits(fn)}, reverse=True)
        classes = sorted({c.strip() for fn in core.function for c in fn.split(";") if c.strip().lower().startswith("bacteriocin")})
        rows.append({"accession": q2acc.get(query, "?"), "contig": query, "AOI": "AOI_" + aoi.split(".")[0],
                     "region_size_bp": int(t.region_size.iloc[0]) if len(t) else 0, "core_peptide_orfs": len(core),
                     "best_hit": hits[0][1] if hits else "motif/HMM only", "all_hits": "; ".join(f"{n} ({sc})" for sc, n in hits),
                     "classes": "; ".join(classes), "session": session})
cols = ["accession", "contig", "AOI", "region_size_bp", "core_peptide_orfs", "best_hit", "all_hits", "classes", "session"]
df = pd.DataFrame(rows, columns=cols)
df = df[df.accession.isin(G)].sort_values(["accession", "AOI"])
df.to_csv(f"{OUT}/bagel4_summary.tsv", sep="\t", index=False)
cm = pd.concat(maps) if maps else pd.DataFrame(columns=["accession"])
cm.to_csv(f"{OUT}/contig_map.tsv", sep="\t", index=False)
screened = set(cm.accession) & set(G)
open(f"{OUT}/coverage.txt", "w").write(f"{len(screened)}/{len(G)} genomes screened; missing: {sorted(set(G) - screened)}\n")
print(f"{len(df)} AOIs; {len(screened)}/{len(G)} genomes screened")
