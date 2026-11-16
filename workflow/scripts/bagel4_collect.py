"""Collect BAGEL4 web-server sessions (bagel4.molgenrug.nl has no API) and summarise areas of interest (AOIs) per genome.
Core peptide = an AOI GeneTable row BAGEL4 colours as core (green) or that carries a scored database hit."""
import argparse, os, urllib.request
from html.parser import HTMLParser
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--sessions", nargs="*", default=[]); ap.add_argument("--genomes", required=True); ap.add_argument("--outdir", required=True)
args = ap.parse_args()
OUT = args.outdir
os.makedirs(f"{OUT}/raw", exist_ok=True)
G = [l.strip() for l in open(args.genomes) if l.strip()]

def db_hits(function):
    """BAGEL4 writes 'score;Name' pairs glued together ('74.1;Mutacin_II166.2;Plantaricin_423').
    Return [(score, name)]: each numeric token scores the next token; a float glued to a name's tail belongs to the next pair."""
    toks = [t.strip() for t in function.split(";")]
    out, score = [], None
    def tail_float(t):
        for i in range(len(t)):
            try:
                v = float(t[i:])
                return (t[:i], v) if "." in t[i:] and i > 0 else None
            except ValueError:
                continue
        return None
    for t in toks:
        if score is not None:
            split = tail_float(t)
            name, nxt = (split if split else (t, None))
            out.append((score, name)); score = nxt; continue
        try:
            score = float(t) if "." in t else None
        except ValueError:
            score = None
    return out


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.files = []
    def handle_starttag(self, tag, attrs):
        href = dict(attrs).get("href", "")
        if tag == "a" and href and not href.startswith(("?", "/")):
            self.files.append(href)




def acc_from_upload(name):
    """BAGEL4 drops the underscore from uploaded file names: GCA061326205.1.fna -> GCA_061326205.1"""
    stem = name.removesuffix(".fna")
    return "GCA_" + stem[3:] if stem.startswith("GCA") and not stem.startswith("GCA_") else stem


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
