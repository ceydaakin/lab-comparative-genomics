"""PHASTEST (phastest.ca) for one genome via its public API, then parse the prophage region table.
The API answers HTTP 500 for multi-FASTA uploads, so contigs are joined with 100-N spacers and regions are mapped back.
The API is intermittently unavailable (empty replies): submission and polling retry with back-off up to --minutes."""
import argparse, json, time, urllib.error, urllib.request
from Bio import SeqIO

API = "https://phastest.ca/phastest_api"
p = argparse.ArgumentParser()
p.add_argument("--fasta", required=True); p.add_argument("--accession", required=True)
p.add_argument("--minutes", type=int, default=240); p.add_argument("--summary", required=True); p.add_argument("--regions", required=True)
a = p.parse_args()
deadline = time.time() + 60 * a.minutes

recs = list(SeqIO.parse(a.fasta, "fasta")); offsets, parts, pos = [], [], 0
for r in recs:
    offsets.append((r.id, pos + 1, pos + len(r.seq))); parts.append(str(r.seq)); pos += len(r.seq) + 100
body = (">joined\n" + ("N" * 100).join(parts) + "\n").encode()


def call(req):
    while True:
        try:
            return json.load(urllib.request.urlopen(req, timeout=180))
        except (urllib.error.URLError, json.JSONDecodeError, ConnectionError, TimeoutError) as e:
            if time.time() > deadline: raise SystemExit(f"PHASTEST unreachable for {a.accession}: {e}")
            time.sleep(120)


job = call(urllib.request.Request(API, data=body, method="POST")).get("job_id")
if not job: raise SystemExit(f"PHASTEST returned no job id for {a.accession}")
while True:
    d = call(urllib.request.Request(f"{API}?acc={job}"))
    if d.get("status") == "Complete": break
    if time.time() > deadline: raise SystemExit(f"PHASTEST job {job} for {a.accession} not finished in {a.minutes} min")
    time.sleep(120)
open(a.summary, "w").write(d["summary"])


def to_contig(region):
    s, e = (int(x) for x in region.split("-"))
    for cid, lo, hi in offsets:
        if lo <= s <= hi: return f"{cid}:{s - lo + 1}-{min(e, hi) - lo + 1}"
    return region


with open(a.regions, "w") as out:
    out.write("accession\tjob\tregion\tlength_kb\tcompleteness\tscore\tposition\tcontig_position\tmost_common_phage\n")
    for line in d["summary"].splitlines():
        t = line.split()
        if len(t) > 12 and t[0].isdigit() and t[1].endswith("Kb") and "(" in t[2]:
            phage = next((x for x in t if x.startswith("PHAGE_")), "").split(",")[0].split("(")[0]
            out.write(f"{a.accession}\t{job}\t{t[0]}\t{t[1][:-2]}\t{t[2].split('(')[0]}\t{t[2].split('(')[1].rstrip(')')}\t{t[4]}\t{to_contig(t[4])}\t{phage}\n")
