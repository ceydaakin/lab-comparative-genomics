"""Translate Roary pan_genome_reference.fa (nucleotide) to protein with the bacterial code (table 11); keeps the header id."""
import sys
from Bio import SeqIO

with open(sys.argv[2], "w") as out:
    for r in SeqIO.parse(sys.argv[1], "fasta"):
        s = r.seq[: len(r.seq) // 3 * 3]
        out.write(f">{r.id}\n{str(s.translate(table=11)).rstrip('*').replace('*', 'X')}\n")
