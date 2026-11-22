# cryptobenthic-gulf-genomics

**From genomes to species assemblages: neutral, selective and functional processes structuring cryptobenthic reef fish diversity in the Gulf of California**

Code, analyses and documentation for project CBF-2026-4124 (SECIHTI, Ciencia Básica y de Frontera 2026), hosted at the Centro de Investigaciones Biológicas del Noroeste (CIBNOR), La Paz, Baja California Sur, Mexico.

> *Español:* Repositorio del proyecto CBF-2026-4124. Contiene código, análisis y documentación sobre la congruencia entre la diversidad de especies y la diversidad genética de peces criptobentónicos arrecifales del Golfo de California, con eDNA, ddRAD, genomas completos y RNA-seq.

## Scientific scope

The project tests whether spatial turnover in community composition and genetic divergence among populations are congruent (species–genetic diversity correlations, SGDCs), and whether that congruence is driven by neutral, selective or functional processes. Cryptobenthic reef fishes of the Gulf of California are the model assemblage.

| Objective | Approach | Folder |
| --- | --- | --- |
| 1. Community structure (alpha and beta diversity) | Anaesthetic collections and eDNA metabarcoding | `analyses/01_community_edna` |
| 1a. eDNA blind spots: reference gaps, primer bias, resolution | In silico analysis of public databases | `analyses/01a_reference_gaps` |
| 1b. Cryptobenthic-specific mitochondrial marker | Mitogenome comparison and primer design | `analyses/01b_marker_design` |
| 2. Neutral and adaptive divergence across species | ddRAD-seq | `analyses/02_ddrad_neutral_adaptive` |
| 3. Reference genomes of three dominant species | PacBio HiFi and Hi-C | `analyses/03_reference_genomes` |
| 4. Genomic regions under selection | High-coverage whole-genome resequencing (HiWGS) | `analyses/04_hiwgs_selection` |
| 5. Function of genes under selection | Annotation, GO and KEGG enrichment | `analyses/05_functional_annotation` |
| 6. Functional expression patterns | RNA-seq | `analyses/06_rnaseq_expression` |
| 7. Hypothesis tests | Mantel tests and distance-decay analyses | `analyses/07_hypothesis_tests` |

## Repository structure

```
analyses/      one folder per objective; each has its own README, code and small result tables
data/          sample metadata only; raw sequence data are deposited in NCBI SRA (see data/README.md)
protocols/     field and laboratory protocols
scripts/       shared utility functions
```

## Data policy

Raw sequencing data are not stored in this repository. They will be deposited in NCBI SRA under a BioProject linked to each publication. Large intermediate files are excluded through `.gitignore`.

## License

Code is released under the MIT License (`LICENSE`). Data, result tables, figures and documentation are released under CC BY 4.0 (`LICENSE-DATA.md`).

## Status

Under active development. Result tables are added to each analysis folder when the corresponding manuscript is submitted.

## Contact

Tania Valdivia Carrillo, CIBNOR. Principal investigator of the project: Fausto Valenzuela Quiñonez, CIBNOR.
