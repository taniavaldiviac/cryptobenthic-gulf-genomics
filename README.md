<p align="center">
  <img src="assets/banner.svg" alt="Cryptobenthic fishes of the Gulf of California: community structure from collections and environmental DNA" width="100%">
</p>

<p align="center">
  <img alt="Status" src="https://img.shields.io/badge/status-in%20development-e07a5f">
  <img alt="Code license" src="https://img.shields.io/badge/code-MIT-3fa7b5">
  <img alt="Data license" src="https://img.shields.io/badge/data-CC%20BY%204.0-3fa7b5">
  <img alt="Python" src="https://img.shields.io/badge/python-3.9%2B-0b3c5d">
  <img alt="Region" src="https://img.shields.io/badge/region-Gulf%20of%20California-0a2a43">
</p>

<p align="center">
  <b>Who lives in the cracks of the reef, and can environmental DNA see them?</b>
</p>

---

Cryptobenthic reef fishes are small, hidden and astonishingly diverse. In the Gulf of California they include gobies, blennies, triplefins, clingfishes and their relatives, many of them endemic. Because they live inside the reef rather than above it, visual censuses miss most of them, and standard fish eDNA markers detect only a handful.

This repository holds the code and analyses for **Objective 1 of project CBF-2026-4124**: describing the spatial structure of the cryptobenthic assemblage across the Gulf of California, and building the molecular tools needed to detect it with environmental DNA.

> **En español.** Repositorio del Objetivo 1 del proyecto CBF-2026-4124: estructura comunitaria espacial (diversidad alfa y beta) de los peces criptobentónicos arrecifales del Golfo de California, a partir de colectas con anestésico y de ADN ambiental, junto con las herramientas moleculares necesarias para detectarlos.

## What this objective does

<table>
  <tr>
    <td width="33%" valign="top">
      <h3>Community structure</h3>
      Alpha and beta diversity of the cryptobenthic assemblage along the Gulf, from anaesthetic collections and eDNA metabarcoding.
      <br><br><code>analyses/community_structure</code>
    </td>
    <td width="33%" valign="top">
      <h3>eDNA blind spots</h3>
      How complete are public reference databases for cryptobenthic species, and how well do standard fish primers amplify and resolve them?
      <br><br><code>analyses/reference_gaps</code>
    </td>
    <td width="33%" valign="top">
      <h3>A marker for the hidden half</h3>
      A cryptobenthic-specific mitochondrial marker designed from new mitogenome references for the Gulf assemblage.
      <br><br><code>analyses/marker_design</code>
    </td>
  </tr>
</table>

## Workflow

```mermaid
flowchart LR
    subgraph Field
        A[Anaesthetic collections<br/>1 m² quadrats] 
        B[Seawater samples<br/>eDNA filtration]
    end
    subgraph Laboratory
        C[Species identification]
        D[Mitogenome sequencing<br/>genome skimming]
        E[Metabarcoding libraries]
    end
    subgraph Analysis
        F[Local reference database]
        G[Taxonomic assignment]
        H[Alpha and beta diversity]
    end
    A --> C --> H
    A --> D --> F
    F --> M[Marker design] --> E
    B --> E --> G
    F --> G --> H
```

## Repository map

```
.
├── analyses/
│   ├── community_structure/   alpha and beta diversity from collections and eDNA
│   ├── reference_gaps/        reference coverage, primer mismatches and resolution (in silico)
│   └── marker_design/         cryptobenthic-specific mitochondrial marker
├── data/                      sample metadata only; sequences go to NCBI
├── protocols/                 field and laboratory protocols
├── scripts/                   shared utilities
└── assets/                    figures for this page
```

## Getting started

```bash
git clone https://github.com/taniavaldiviac/cryptobenthic-gulf-genomics.git
cd cryptobenthic-gulf-genomics
pip install -r requirements.txt

# reference coverage, primer and resolution analysis
python analyses/reference_gaps/cripticos_insilico.py --email you@example.org --outdir analyses/reference_gaps/results
```

An NCBI API key (`--api-key`) speeds up queries from 3 to 10 requests per second.

## Data

Raw sequence data are not stored here. eDNA reads and new mitogenome references will be deposited in NCBI under a BioProject linked to each publication; accessions will be listed in [`data/README.md`](data/README.md). Result tables are added to each analysis folder when the corresponding manuscript is submitted.

## Project

<table>
  <tr>
    <td valign="top"><b>Project</b></td>
    <td>Del genoma al ensamblaje de especies: evaluación de los procesos evolutivos neutrales, selectivos y funcionales que estructuran la biodiversidad en ambientes marinos (CBF-2026-4124)</td>
  </tr>
  <tr>
    <td valign="top"><b>Funding</b></td>
    <td>SECIHTI, Ciencia Básica y de Frontera 2026</td>
  </tr>
  <tr>
    <td valign="top"><b>Host institution</b></td>
    <td>Centro de Investigaciones Biológicas del Noroeste (CIBNOR), La Paz, Baja California Sur, Mexico</td>
  </tr>
  <tr>
    <td valign="top"><b>Objective 1 lead</b></td>
    <td>Tania Valdivia Carrillo</td>
  </tr>
  <tr>
    <td valign="top"><b>Principal investigator</b></td>
    <td>Fausto Valenzuela Quiñonez</td>
  </tr>
</table>

## License

Code: [MIT](LICENSE). Data, result tables, figures and documentation: [CC BY 4.0](LICENSE-DATA.md).
