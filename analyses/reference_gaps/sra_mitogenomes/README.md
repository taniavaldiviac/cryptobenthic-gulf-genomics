# Mitogenomes recovered from public raw reads (SRA/ENA)

Recovers mitochondrial reference sequences for Gulf of California fishes that lack them, using raw reads already
deposited in the Sequence Read Archive. The goal is to fill gaps in the 12S, 16S and COI regions targeted by common
eDNA metabarcoding primers (MiFish-U, teleo, Riaz 12S-V5, MarVer1, Fish16S, Leray COI).

## Workflow
1. **Search.** We query SRA for each species of the regional checklist, using the accepted name and the GBIF synonym.
   - **Retained:** WGS and targeted-capture runs (UCE, exon capture).
   - **Exploratory only:** ddRAD.
   - **Discarded:** amplicon records.
2. **Download.** Paired reads come from ENA.
   - **WGS:** the first 5 × 10⁶ reads of each file. This is enough for the high-copy mitochondrial genome.
   - **Capture:** the full run.
   - **Check:** if R1 and R2 arrive with different read counts (truncated transfers), both are trimmed to the shared length.
3. **Assembly.** GetOrganelle v1.7.7.1 (`-F animal_mt`).
   - **Seed:** the closest available mitogenome (same genus or family), or a multi-teleost seed.
   - **Assembler:** SPAdes v4.3.0 in `--only-assembler` mode (k = 21, 55, 85), with memory limits for small machines.
4. **Fallback for low mitochondrial coverage.** BWA-MEM v0.7.19 mapping to the seed, then a samtools v1.25 consensus (minimum depth 2–3×).
5. **Validation and marker extraction.**
   - **Identity:** BLAST against conspecific GenBank fragments.
   - **Markers:** `eval_markers.py` locates each primer pair, IUPAC-aware with ≤4 mismatches, and treats complete sequences as circular.
   - **Primer sites:** `primer_mismatch_detail.py` prints mismatch positions at each site.

## Files
| File | Purpose |
|---|---|
| `sra_runs.tsv` | Runs evaluated, with strategy, assembly mode and seed |
| `run_one.sh` | Download from ENA + GetOrganelle or mapping for one run |
| `run_batch.sh` | Sequential driver over a TSV (species, run, seed, strategy) |
| `sra_mitogenomes_local.sh` | Stand-alone version for a workstation or cluster (sra-tools) |
| `eval_markers.py` | Extracts marker amplicons and primer mismatches from an assembly |
| `primer_mismatch_detail.py` | Per-site mismatch map for each primer |
| `environment.yml` | Conda environment |

```bash
conda env create -f environment.yml && conda activate sramito
get_organelle_config.py -a animal_mt
bash run_batch.sh batch.tsv
python eval_markers.py asm/<species_run>/*path_sequence.fasta amplicons.fa
```

## Notes
- **Results:** assemblies and summary tables are kept in `results/`, which is not tracked. They will be released upon manuscript submission.
- **Reference-bias caveats:**
  - Consensus sequences from mapping, and assemblies longer than expected, are treated as provisional.
  - ddRAD libraries gave no usable coverage of 12S, 16S or COI.
- **Data attribution:** reads belong to their original BioProjects, which must be cited (see `sra_runs.tsv` and the run accessions).
