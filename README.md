# Protein Similarity & Functional Conservation Analysis

An end-to-end computational biology and comparative bioinformatics pipeline for analyzing evolutionary divergence, sequence homology, information-theoretic conservation, and active site resilience across protein families.

![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Biopython](https://img.shields.io/badge/Biopython-1.88-orange.svg)
![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)

---

## Overview

In comparative proteomics and protein engineering, high sequence similarity does not always imply identical function, nor does low sequence identity (the evolutionary "twilight zone", ~20–35% identity) necessarily indicate functional divergence. Critical catalytic residues often remain strictly invariant under intense selective pressure even when background sequences drift substantially.

This toolkit provides a complete, modular pipeline bridging **primary sequence alignment**, **information theory (Shannon entropy)**, **phylogenetic inference**, and **3D macromolecular visualization**:

```
                         ┌──────────────────────────────────────────────┐
                         │           Input FASTA Sequences              │
                         └──────────────────────┬───────────────────────┘
                                                │
                     ┌──────────────────────────┴──────────────────────────┐
                     ▼                                                     ▼
    ┌─────────────────────────────────┐                   ┌─────────────────────────────────┐
    │       Adaptive Pairwise         │                   │    Progressive Multiple Seq     │
    │      Alignment & Distance       │                   │         Alignment (MSA)         │
    │  (Needleman-Wunsch / S-Waterman)│                   └────────────────┬────────────────┘
    └────────────────┬────────────────┘                                    │
                     ▼                                                     ▼
    ┌─────────────────────────────────┐                   ┌─────────────────────────────────┐
    │   Hierarchical Clustering &     │                   │     Shannon Entropy Column      │
    │  Phylogenetic Trees (NJ/UPGMA)  │                   │       Conservation (0.0 - 1.0)  │
    └─────────────────────────────────┘                   └────────────────┬────────────────┘
                                                                           │
                     ┌─────────────────────────────────────────────────────┘
                     ▼
    ┌───────────────────────────────────────────────────────────────────┐
    │          Active Site & Functional Tolerance Analysis              │
    │   - Universal Coordinate Mapping (1-indexed reference)            │
    │   - Physicochemical Divergence Checks (Charge / Polarity / Size)  │
    │   - 3D Structural Inspection & Py3Dmol Visualization (PDB)        │
    └───────────────────────────────────────────────────────────────────┘
```

---

## Key Features

1. **Deterministic Adaptive Alignment Engine** (`src/aligner.py`):
   - **Dynamic Algorithm Routing**: Automatically switches between **Needleman-Wunsch** (global) and **Smith-Waterman** (local) based on pairwise sequence length asymmetry:
     $$\text{Length Ratio} = \frac{\min(L_1, L_2)}{\max(L_1, L_2)}$$
     - If $\text{Ratio} \ge 0.75 \implies$ Global alignment.
     - If $\text{Ratio} < 0.75 \implies$ Local alignment (avoids distortion in multi-domain proteins with flanking variable extensions).
   - Standard BLOSUM62 scoring matrix with Gotoh affine gap penalties (gap open: $-10.0$, gap extend: $-0.5$).
   - Symmetric pairwise percent identity and evolutionary distance matrices ($100 - \% \text{ID}$).

2. **Progressive Multiple Sequence Alignment (MSA)** (`src/msa.py`):
   - Constructs uniform profile alignments across diverse homologous families while strictly preserving residue integrity and eliminating alignment drift.

3. **Shannon Entropy & Conservation Profiling** (`src/entropy.py`):
   - Computes per-column information-theoretic Shannon entropy ($H$ in bits):
     $$H = -\sum_{i=1}^{20} p_i \log_2(p_i)$$
   - Calculates normalized conservation index ($1 - H / H_{\max}$, where $H_{\max} = \log_2(20) \approx 4.32$ bits).
   - Generates publication-grade conservation plots highlighting catalytic columns.

4. **Universal Active Site & Physicochemical Tolerance Mapping** (`src/active_site_analyzer.py`):
   - Ingests arbitrary JSON active site specifications mapped to a 1-indexed reference protein.
   - Evaluates mutations against amino acid physicochemical properties:
     - **Identical**: Exact residue conservation.
     - **Conservative Swap**: Same charge/polarity group (e.g. $\text{Glu} \to \text{Asp}$, $\text{Lys} \to \text{Arg}$).
     - **Diverged**: Radical shift in charge, polarity, or hydrophobicity (e.g. $\text{Acidic} \to \text{Aliphatic}$).

5. **Phylogenetic Tree Reconstruction** (`src/phylogeny.py`):
   - Reconstructs evolutionary trees using **Neighbor-Joining (NJ)** and **UPGMA** from pairwise distance matrices.
   - Computes explicit branch lengths and exports standard Newick (`.nwk`) formats.

6. **UniProt REST API Integration** (`src/annotation_fetcher.py`):
   - Automatically downloads curated, experimentally validated active sites, binding sites, and metal-binding residues from UniProtKB REST endpoints.

7. **3D Structural Mapping & Interactive Visualization** (`src/structure_viewer.py`):
   - Directly downloads PDB coordinates from RCSB.
   - Renders interactive 3D cartoon backbones with highlighted catalytic residues in stick representations via `py3Dmol`.

---

## Project Structure

```text
Protein Similarity Analysis/
├── config/
│   └── sites/                    # Functional residue specifications
│       ├── globins.json          # Distal and Proximal Histidines (HBA_HUMAN)
│       ├── mpro.json             # SARS-CoV-2 Main Protease catalytic dyad
│       ├── uniprot_trypsin.json  # Serine protease catalytic triad & binding pocket
│       └── uniprot_pka.json      # Kinase catalytic core, DFG motif & salt bridge
├── data/
│   └── raw/                      # Curated benchmark FASTA datasets
│       ├── sequences.fasta       # Globin superfamily (HBA, HBB, MYG, NGB)
│       ├── proteases.fasta       # Serine proteases (Trypsin, Chymotrypsin, Elastase, etc.)
│       └── kinases.fasta         # Protein Kinases (PKA, Src, EGFR, CDK2, etc.)
├── notebooks/                    # In-depth biological case studies
│   ├── 01_globin_conservation.ipynb
│   ├── 02_serine_protease_triad.ipynb
│   └── 03_kinase_catalytic_core.ipynb
├── results/                      # Generated figures, tables, and Newick trees
│   ├── figures/                  # Heatmaps, dendrograms, entropy profiles
│   ├── tables/                   # Identity matrices, conservation metrics, site summaries
│   └── trees/                    # Exported Newick (.nwk) phylogenetic trees
├── src/                          # Core analytical modules
│   ├── __init__.py
│   ├── active_site_analyzer.py   # Universal active site mapper & tolerance classifier
│   ├── aligner.py                # Adaptive Needleman-Wunsch / Smith-Waterman aligner
│   ├── annotation_fetcher.py     # UniProtKB REST API feature parser
│   ├── clustering.py             # UPGMA linkage & hierarchical clustering
│   ├── data_loader.py            # FASTA parser & IUPAC sequence validation
│   ├── entropy.py                # Shannon entropy & conservation scoring
│   ├── msa.py                    # Progressive multiple sequence alignment
│   ├── phylogeny.py              # Neighbor-Joining & UPGMA tree reconstruction
│   ├── structure_viewer.py       # 3D Py3Dmol PDB macromolecular renderer
│   └── visualize.py              # Publication-grade plotting utilities
├── tests/                        # Automated unit tests
│   ├── __init__.py
│   ├── test_active_analyzer.py   # Active site & divergence tests
│   └── test_aligner.py           # Alignment & matrix unit tests
├── main.py                       # CLI orchestrator for end-to-end analysis
├── run_active_sites.py           # CLI tool for site conservation analysis
├── requirements.txt              # Pinned Python package dependencies
└── README.md
```

---

## Installation & Setup

### 1. Clone the Repository
```bash
git clone https://github.com/AmaLiVeAmaLiVe/Protein-Similarity-Analysis.git
cd "Protein Similarity Analysis"
```

### 2. Set Up a Virtual Environment

**On Windows (PowerShell):**
```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```
*(If script execution is disabled, run `Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned` first).*

**On macOS / Linux:**
```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

---

## Usage

### 1. Run the Full End-to-End Pipeline via CLI
Execute the automated 5-stage pipeline on the benchmark globin dataset:

```bash
python main.py --input data/raw/sequences.fasta --sites-config config/sites/globins.json --output-dir results
```

**Pipeline Stages:**
```text
[1/5] Loading sequences from: data\raw\sequences.fasta (IUPAC validation)
[2/5] Computing pairwise alignments & distance tree (Adaptive BLOSUM62)
[3/5] Constructing Progressive Multiple Sequence Alignment (MSA)
[4/5] Computing per-position Shannon Entropy profile (Conservation indexing)
[5/5] Performing active site tolerance analysis (Physicochemical classification)

All artifacts generated successfully in: results/
```

### 2. Inspect Custom Active Site Configurations
Analyze any protein family against an external JSON site specification:

```bash
python run_active_sites.py \
    --input data/raw/proteases.fasta \
    --config config/sites/uniprot_trypsin.json \
    --output results/tables/trypsin_conservation.csv
```

### 3. Fetch UniProt Annotations Dynamically
```python
from src.annotation_fetcher import fetch_uniprot_features

# Fetch curated functional sites for Human Hemoglobin Alpha
features = fetch_uniprot_features("HBA_HUMAN", output_path="config/sites/hba_uniprot.json")
```

---

## Biological Case Studies (Notebooks)

The repository includes three comprehensive, runnable case studies in `notebooks/`:

### Case 1: Globin Superfamily (`01_globin_conservation.ipynb`)
- **Proteins**: Human Hemoglobin $\alpha$ (`HBA_HUMAN`), Hemoglobin $\beta$ (`HBB_HUMAN`), Human Myoglobin (`MYG_HUMAN`), Sperm Whale Myoglobin (`MYG_PHYMC`), and Neuroglobin (`NGB_HUMAN`).
- **Biological Finding**: Neuroglobin diverges to the sequence "twilight zone" (~$21-25\%$ identity to hemoglobins/myoglobins). Despite extensive sequence drift, the **Distal Histidine** ($\text{His58/59}$) and **Proximal Histidine** ($\text{His87/88}$) remain strictly conserved, confirming structural and functional heme-binding invariance.

### Case 2: Serine Protease Catalytic Triad (`02_serine_protease_triad.ipynb`)
- **Proteins**: Bovine Trypsin, Human Chymotrypsin, Elastase, Thrombin, Coagulation Factor X, and Plasmin.
- **Biological Finding**: Demonstrates the conservation of the spatial catalytic triad ($\text{His57}-\text{Asp102}-\text{Ser195}$) across enzymes with diverging primary substrate binding pockets ($\text{S1}$ specificity pocket).

### Case 3: Protein Kinase Catalytic Core (`03_kinase_catalytic_core.ipynb`)
- **Proteins**: PKA, Src, EGFR, CDK2, and MAP kinase.
- **Biological Finding**: Demonstrates dynamic local/global alignment routing on multi-domain kinases with large regulatory insertions. Tracks conservation of the **DFG motif** ($\text{Asp}-\text{Phe}-\text{Gly}$), catalytic base ($\text{Asp166}$), and the activating $\text{Lys}-\text{Glu}$ salt bridge.

---

## Testing

Run unit tests using `pytest` within your virtual environment:

```powershell
python -m pytest tests/
```

Test suite covers:
- Verification of percentage identity calculation on identical and mismatched pairs.
- Invariance and symmetry of distance and identity matrices.
- Correct classification of conservative amino acid swaps versus radical divergence.

---

## Scientific References

- **Needleman, S. B., & Wunsch, C. D. (1970).** A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology*, 48(3), 443-453.
- **Smith, T. F., & Waterman, M. S. (1981).** Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195-197.
- **Henikoff, S., & Henikoff, J. G. (1992).** Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences*, 89(22), 10915-10919.
- **Shannon, C. E. (1948).** A mathematical theory of communication. *Bell System Technical Journal*, 27(3), 379-423.
- **Saitou, N., & Nei, M. (1987).** The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution*, 4(4), 406-425.
- **Cock, P. J., et al. (2009).** Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.

---

## License

This project is open-source under the [MIT License](LICENSE).

