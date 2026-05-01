# Project File Structure

```
nhej-repair-genome-analysis/
│
├── README.md
│   └── Project overview and documentation
│
├── data/
│   ├── raw/
│   │   └── Raw bacterial genomes from GTDB and UniProt
│   │
│   └── processed/
│       └── Pre-processed and cleaned datasets
│
├── scripts/
│   │
│   ├── hmm/
│   │   ├── hmm_detection.py
│   │   └── profile_builder.py
│   │
│   ├── interpro/
│   │   ├── annotation.py
│   │   └── result_parser.py
│   │
│   ├── phylogeny/
│   │   ├── mapping.py
│   │   └── tree_analysis.py
│   │
│   └── gc_analysis/
│       ├── gc_content.py
│       └── correlation_analysis.py
│
├── results/
│   │
│   ├── hmm_profiles/
│   │   └── HMM profile definitions and outputs
│   │
│   ├── hmmer/
│   │   └── HMMER tool result files
│   │
│   ├── interpro/
│   │   └── InterPro domain detection results
│   │
│   ├── phylogeny/
│   │   ├── trees/
│   │   │   └── Phylogenetic tree files
│   │   │
│   │   └── mappings/
│   │       └── Phylogenetic mapping results
│   │
│   └── gc/
│       ├── distribution/
│       │   └── GC content distribution analysis
│       │
│       └── correlation/
│           └── Correlation with NHEJ repair systems
│
├── figures/
│   │
│   ├── hmm_results/
│   │   └── HMM-based visualizations
│   │
│   ├── interpro_results/
│   │   └── Domain annotation plots
│   │
│   ├── phylogeny_results/
│   │   └── Phylogenetic tree visualizations
│   │
│   └── gc_analysis/
│       └── GC content distribution plots
│
├── notebooks/
│   │
│   ├── 01_data_exploration.ipynb
│   │   └── Initial data exploration and summary statistics
│   │
│   ├── 02_hmm_analysis.ipynb
│   │   └── HMM detection workflow and results
│   │
│   ├── 03_interpro_analysis.ipynb
│   │   └── InterPro annotation exploration
│   │
│   ├── 04_phylogeny_analysis.ipynb
│   │   └── Phylogenetic mapping and analysis
│   │
│   └── 05_gc_correlation.ipynb
│       └── GC content correlation analysis
│
└── docs/
    │
    ├── methods.md
    │   └── Detailed methodology documentation
    │
    ├── data_sources.md
    │   └── Data source information and citations
    │
    ├── installation.md
    │   └── Installation and setup instructions
    │
    └── results_interpretation.md
        └── Interpretation guide for results
```

## Directory Descriptions

### `data/`
Contains raw and processed datasets used in the analysis pipeline:
- **raw/**: Original bacterial genomes from GTDB and UniProt
- **processed/**: Pre-processed and cleaned datasets ready for analysis

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
- **phylogeny/**: Phylogenetic tree files and mappings
- **gc/**: GC content distribution and correlation analysis

### `figures/`
Final publication-ready plots and visualizations organized by analysis type:
- **hmm_results/**: HMM-based detection visualizations
- **interpro_results/**: Domain annotation plots
- **phylogeny_results/**: Phylogenetic tree visualizations
- **gc_analysis/**: GC content distribution and correlation plots

### `notebooks/`
Jupyter notebooks for exploratory analysis and visualization:
- `01_data_exploration.ipynb`: Initial data exploration and summary statistics
- `02_hmm_analysis.ipynb`: HMM detection workflow and results
- `03_interpro_analysis.ipynb`: InterPro annotation exploration
- `04_phylogeny_analysis.ipynb`: Phylogenetic mapping and analysis
- `05_gc_correlation.ipynb`: GC content correlation analysis

### `docs/`
Supporting documentation:
- **methods.md**: Detailed methodology documentation
- **data_sources.md**: Data source information and citations
- **installation.md**: Installation and setup instructions
- **results_interpretation.md**: Interpretation guide for results

### `README.md`
Project overview, objectives, methodology, and reproducibility information.
