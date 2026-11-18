# Tools, thresholds and their sources

The group's papers:
- **P22:** Gumustop & Ortakcı 2022, *BMC Genomics*, doi:10.1186/s12864-022-09053-y (*L. parabuchneri*).
- **P23:** Gumustop & Ortakcı 2023, *CSBJ*, doi:10.1016/j.csbj.2023.10.004 (*L. coryniformis*).
- **P25:** Gumustop, Genel, Kurt & Ortakcı 2025, *PLoS One*, doi:10.1371/journal.pone.0325832 (*L. buchneri*).
- **S26:** Kasapoglu et al. 2026, *Sci Rep*, doi:10.1038/s41598-026-51219-1 (*L. fermentum*).
- **K26:** Kurt et al. 2026, *BMC Genomics*, doi:10.1186/s12864-026-13098-8 (*L. curvatus*).

| Step | Group method | This workflow | Setting | Note |
|---|---|---|---|---|
| QC | BUSCO, lactobacillales_odb10 | BUSCO 6.1 | complete ≥ 98% (P25, K26); 95% in P22 | `qc.min_busco_complete` |
| ANI | GET_HOMOLOGUES CDS-ANI, 95–96% | skani | species 95%; clonal 99.9% | GET_HOMOLOGUES not installable on current bioconda |
| dDDH | TYGS | TYGS web, collected by job GUID | d4, species cut-off 70% | submission form needs an e-mail |
| Annotation | Prokka 1.14.x | Prokka 1.14.6 | `--kingdom Bacteria` | |
| Pangenome | Roary 3.13.0 (-e -n -v -r), 95% | Roary 3.13.0 | `-i 95 -cd 100` | `-r` plots replaced by our figures |
| Openness | micropan Heaps' law, 10,000 permutations | micropan 2.1 (CRAN) | α < 1 = open | |
| Structure | prabclus Jaccard + PCoA | prabclus + cmdscale | | |
| Phylogeny | Roary core + FastTree GTR (P22, P23); snippy + Gubbins + FastTree (P25, K26) | Roary core + FastTree GTR; Gubbins optional | | snippy not used: the core alignment comes from Roary |
| COG | eggNOG-mapper v2 / PSI-BLAST vs COG (P22) | rpsblast vs NCBI CDD COG, COG 2024 | E < 1e-5 | eggNOG web needs an e-mail |
| CAZymes | dbCAN2, HMMER | dbCAN HMM db v5.2.9 + pyhmmer | E < 1e-15, coverage > 0.35 (P22) | |
| Bacteriocins | BAGEL4 | BAGEL4 web, collected by session id | | no API |
| Prophages | PHASTER / PHASTEST | PHASTEST API | intact > 90, questionable 70–90 | contigs joined with 100 N |
| CRISPR | CRISPRCasFinder, CRISPRviz, CRISPRClassify | MinCED + Cas typing (Makarova 2020 signatures) | MMseqs2 ≥ 50% id, ≥ 70% cov to PGAP Cas | Prokka names miss type I-B Cas proteins |
| AMR | CARD strict/perfect | RGI 6 + CARD | Strict/Perfect | |
| Plasmids | PLSDB web screen | mash screen -w vs PLSDB 2024_05_31_v2 | identity ≥ 0.99, p ≤ 0.1; present at ≥ 50% hashes | |
| IS elements | ISfinder | ISEScan | single-threaded | multi-threaded FragGeneScan fails silently on macOS |
| Comparison | BRIG, 70/90% rings | BLASTn rings, BRIG method | 1 kb bins | BRIG is a desktop GUI |
| Genomic islands | GIPSy | focal-only regions from BRIG BLASTn (≥ 5 kb, < 70% in ≥ 90% genomes) | | GIPSy is a desktop GUI |
