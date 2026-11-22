# cryptobenthic-gulf-genomics

**Community structure of cryptobenthic reef fishes in the Gulf of California from anaesthetic collections and environmental DNA**

Code, analyses and documentation for Objective 1 of project CBF-2026-4124 (SECIHTI, Ciencia Básica y de Frontera 2026), hosted at the Centro de Investigaciones Biológicas del Noroeste (CIBNOR), La Paz, Baja California Sur, Mexico.

> *Español:* Repositorio del Objetivo 1 del proyecto CBF-2026-4124: estructura comunitaria espacial (diversidad alfa y beta) de los peces criptobentónicos arrecifales del Golfo de California a partir de colectas con anestésico y de ADN ambiental (eDNA).

## Scope

Objective 1 characterizes the spatial community structure (alpha and beta diversity) of the cryptobenthic reef fish assemblage across the Gulf of California, combining anaesthetic collections with eDNA metabarcoding. Because standard fish metabarcoding markers detect few cryptobenthic species, this objective also evaluates public reference databases and primers, and develops a cryptobenthic-specific mitochondrial marker.

| Component | Approach | Folder |
| --- | --- | --- |
| Community structure (alpha and beta diversity) | Anaesthetic collections and eDNA metabarcoding | `analyses/community_structure` |
| Reference databases, primer bias and taxonomic resolution | In silico analysis of public databases | `analyses/reference_gaps` |
| Cryptobenthic-specific mitochondrial marker | Mitogenome comparison and primer design | `analyses/marker_design` |

## Repository structure

```
analyses/      one folder per component; each has its own README, code and small result tables
data/          sample metadata only; raw sequence data are deposited in NCBI (see data/README.md)
protocols/     field and laboratory protocols
scripts/       shared utility functions
```

## Data policy

Raw sequencing data are not stored in this repository. They will be deposited in NCBI under a BioProject linked to each publication. Large intermediate files are excluded through `.gitignore`.

## License

Code is released under the MIT License (`LICENSE`). Data, result tables, figures and documentation are released under CC BY 4.0 (`LICENSE-DATA.md`).

## Status

Under active development. Result tables are added to each analysis folder when the corresponding manuscript is submitted.

## Contact

Tania Valdivia Carrillo, CIBNOR. Principal investigator of the project: Fausto Valenzuela Quiñonez, CIBNOR.
