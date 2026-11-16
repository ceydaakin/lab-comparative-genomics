# lab-comparative-genomics

A reproducible Snakemake workflow for population-wide comparative genomics of lactic acid bacteria. It encodes the pipeline used in the Ortakcı group's papers, so a new species needs only a new genome list and a few lines of config. An optional step has Claude draft the Results text from the computed tables.

First applied to *Liquorilactobacillus nagelii* (46 NCBI genomes, focal strain AGA58 from shalgam).

## What it runs

| Step | Tool (as in the group's papers) | Output |
|---|---|---|
| Download | NCBI Datasets CLI | genomes, PGAP proteins, assembly metadata |
| QC | BUSCO 6, lineage and cut-off from config (≥98%) | `results/qc/`, genome set (checkpoint) |
| Species, clones | skani ANI; single-linkage clonal groups | `results/ani/` |
| Annotation | Prokka 1.14.6 (`--kingdom Bacteria`) | `results/prokka/` |
| Pangenome | Roary 3.13.0 (95% identity, core = 100%); micropan Heaps' law; prabclus Jaccard PCoA | `results/roary/`, `results/pangenome/` |
| Phylogeny | Roary core alignment → FastTree 2 GTR; Gubbins recombination filter (optional) | `results/phylogeny/` (iTOL-ready Newick) |
| COG | rpsblast vs NCBI CDD COG, COG 2024 categories | `results/cog/` |
| CAZymes | dbCAN HMMs (E < 1e-15, coverage > 0.35) | `results/cazy/` |
| CRISPR-Cas | MinCED arrays; Cas typing against NCBI PGAP Cas proteins | `results/crispr/` |
| Prophages | PHASTEST (public API) | `results/phastest/` |
| Bacteriocins | BAGEL4 (web; sessions collected by id) | `results/bagel4/` |
| IS elements | ISEScan (ISfinder families) | `results/isescan/` |
| AMR | RGI + CARD (Strict/Perfect) | `results/rgi/` |
| Plasmids | mash screen vs PLSDB 2024_05_31_v2 | `results/plsdb/` |
| Genome comparison | BLASTn rings in the BRIG manner against the focal strain | `results/figures/fig10_brig.*`, `results/brig/` |
| Report | Figures (PNG 300 dpi + PDF), `Supplementary_Tables.xlsx`, `key_numbers.json` | `results/figures/`, `results/report/` |
| AI draft (optional) | Claude writes Results + figure legends from `key_numbers.json` only | `results/report/results_draft.md` |

`docs/tools.md` maps each step to the paper it comes from and gives every threshold.

## Run

Requirements: Linux (or WSL2) with conda/mamba, and Snakemake 8 or 9. About 20 GB of disk for the databases and conda environments.

```bash
mamba create -n snakemake -c conda-forge -c bioconda "snakemake>=8,<10"
mamba activate snakemake
snakemake --sdm conda --cores 16 -n        # dry run: shows every job
snakemake --sdm conda --cores 16 --resources phastest_jobs=4
```

To run a new species:
1. Put the NCBI assembly accessions in `config/genomes.tsv`.
2. In `config/config.yaml`, set `focal`, `species`, the BUSCO lineage and the `niches` keywords.
3. Run Snakemake as above.

Test run (5 genomes, PHASTEST off): `snakemake --sdm conda --cores 8 --configfile .tests/config_test.yaml`.

### Steps that use web services
- **PHASTEST:** genomes are submitted to the public API at phastest.ca. Concurrency is capped with `--resources phastest_jobs=N`; switch it off with `phastest.enabled: false`.
- **BAGEL4:** has no API.
  1. Upload `resources/genomes/*.fna` on bagel4.molgenrug.nl (40 genomes or fewer per session).
  2. Copy each session id (the page's `sessionID` variable) into `bagel4.sessions`.
  3. Rerun; the pipeline collects and parses the results.

### Optional AI report
Set `ai_report.enabled: true`. Credentials come from `ANTHROPIC_API_KEY` or an `ant auth login` profile. Claude receives only `key_numbers.json` (numbers the pipeline computed) and is told not to cite anything else. The output is a draft for expert review, marked as such. Each run is a billed API call. `python workflow/scripts/ai_report.py ... --dry-run` prints the prompt without calling the API.

## Known limitations
- **Apple Silicon:** Prokka and Gubbins have no working builds there, so run on Linux. Docker works on a Mac via `--platform linux/amd64`.
- **Tools not reproduced locally:** TYGS (dDDH), eggNOG-mapper/KAAS and CRISPRCasFinder are web forms that need an e-mail address or have no API. GET_HOMOLOGUES does not install on current bioconda; skani ANI is used instead.
- **PHASTEST input:** the API rejects multi-FASTA uploads, so contigs are joined with 100-N spacers and region coordinates are mapped back to the contigs.

## Layout
Follows the Snakemake deployment guidelines: `workflow/` (Snakefile, `rules/`, `envs/`, `scripts/`), `config/`, `results/`, `resources/`.
