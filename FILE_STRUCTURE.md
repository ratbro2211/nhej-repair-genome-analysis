# Project File Structure

```
nhej-repair-genome-analysis/
│
├── data/
│   └── Raw and processed datasets
│
├── scripts/
│   ├── hmm/
│   │   └── HMM-based detection scripts
│   │
│   ├── interpro/
│   │   └── InterProScan annotation scripts
│   │
│   ├── phylogeny/
│   │   └── Phylogenetic analysis scripts
│   │
│   └── gc_analysis/
│       └── GC content analysis scripts
│
├── results/
│   ├── hmm_profiles/
│   │   └── HMM profile results
│   │
│   ├── hmmer/
│   │   └── HMMER output files
│   │
│   ├── interpro/
│   │   └── InterPro annotation results
│   │
│   ├── phylogeny/
│   │   └── Phylogenetic analysis results
│   │
│   └── gc/
│       └── GC content analysis results
│
├── figures/
│   └── Final plots and visualizations
│
├── notebooks/
│   └── Jupyter notebooks (if any)
│
├── docs/
│   └── Supporting documentation
│
└── README.md
    └── Project overview and documentation
```

## Directory Descriptions

### `data/`
Contains raw and processed datasets used in the analysis pipeline, including bacterial genomes from GTDB and UniProt.

### `scripts/`
Main analysis scripts organized by methodology:
- **hmm/**: Hidden Markov Model-based detection scripts
- **interpro/**: InterProScan annotation and processing scripts
- **phylogeny/**: Phylogenetic mapping and analysis scripts
- **gc_analysis/**: GC content comparison and analysis scripts

### `results/`
Output files and analysis results organized by method:
- **hmm_profiles/**: HMM profile definitions and outputs
- **hmmer/**: HMMER tool results
- **interpro/**: InterPro domain detection results
- **phylogeny/**: Phylogenetic tree and mapping results
- **gc/**: GC content distribution analysis results

### `figures/`
Final publication-ready plots and visualizations.

### `notebooks/`
Jupyter notebooks for exploratory analysis and visualization.

### `docs/`
Supporting documentation, methods descriptions, and additional references.

### `README.md`
Project overview, objectives, methodology, and reproducibility information.
