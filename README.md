# nhej-repair-genome-analysis
Large-scale analysis of NHEJ repair systems and their association with GC content across bacterial genomes.

# 🧬 GC Content and Evolution of NHEJ Repair Systems in Bacteria

## Overview
This project investigates the distribution and evolution of bacterial Non-Homologous End Joining (NHEJ) DNA repair systems and their association with genomic GC content. The study integrates large-scale genomic datasets with phylogenetic analysis to understand the evolutionary dynamics of DNA repair mechanisms across bacterial diversity.

---

## Research Objective
The primary goal of this study is to understand the biological and evolutionary factors contributing to GC content diversity in bacterial genomes, with a particular focus on the role of DNA repair systems, particularly NHEJ components (Ku and LigD), in shaping genomic composition.

---

## Key Questions
- Is the presence of NHEJ components (Ku and LigD) associated with GC-rich genomes?
- How are DNA repair systems distributed across bacterial phylogeny?
- Can large-scale datasets provide more robust insights compared to previous studies?

---

## Methodology Overview

The analysis pipeline consists of the following steps:

1. **Data Collection**
   - Bacterial genomes obtained from GTDB and UniProt reference proteomes

2. **Domain Detection**
   - HMM-based detection (high sensitivity)
   - InterProScan annotation (high specificity)

3. **Filtering and Processing**
   - Removal of incomplete or low-confidence hits
   - Selection of common genome sets for comparison

4. **Method Comparison**
   - Comparison between HMM results, InterPro annotations, and published datasets

5. **Phylogenetic Mapping**
   - Mapping of Ku presence/absence onto bacterial phylogenetic trees

6. **GC Content Analysis**
   - Comparative analysis of GC content distributions
   - Visualization using violin plots

---

## Repository Structure

```
nhej-repair-genome-analysis/
│
├── README.md
│   └── Project overview and documentation
│
├── data/
│   │
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
│   │   │
│   │   ├── trees/
│   │   │   └── Phylogenetic tree files
│   │   │
│   │   └── mappings/
│   │       └── Phylogenetic mapping results
│   │
│   └── gc/
│       │
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

### Directory Descriptions

**data/**
- Contains raw and processed datasets used in the analysis pipeline
- `raw/`: Original bacterial genomes from GTDB and UniProt
- `processed/`: Pre-processed and cleaned datasets ready for analysis

**scripts/**
- Main analysis scripts organized by methodology
- `hmm/`: Hidden Markov Model-based detection scripts
- `interpro/`: InterProScan annotation and processing scripts
- `phylogeny/`: Phylogenetic mapping and analysis scripts
- `gc_analysis/`: GC content comparison and analysis scripts

**results/**
- Output files and analysis results organized by method
- `hmm_profiles/`: HMM profile definitions and outputs
- `hmmer/`: HMMER tool results
- `interpro/`: InterPro domain detection results
- `phylogeny/`: Phylogenetic tree files and mappings
- `gc/`: GC content distribution and correlation analysis

**figures/**
- Final publication-ready plots and visualizations organized by analysis type
- `hmm_results/`: HMM-based detection visualizations
- `interpro_results/`: Domain annotation plots
- `phylogeny_results/`: Phylogenetic tree visualizations
- `gc_analysis/`: GC content distribution and correlation plots

**notebooks/**
- Jupyter notebooks for exploratory analysis and visualization
- `01_data_exploration.ipynb`: Initial data exploration and summary statistics
- `02_hmm_analysis.ipynb`: HMM detection workflow and results
- `03_interpro_analysis.ipynb`: InterPro annotation exploration
- `04_phylogeny_analysis.ipynb`: Phylogenetic mapping and analysis
- `05_gc_correlation.ipynb`: GC content correlation analysis

**docs/**
- Supporting documentation
- `methods.md`: Detailed methodology documentation
- `data_sources.md`: Data source information and citations
- `installation.md`: Installation and setup instructions
- `results_interpretation.md`: Interpretation guide for results

---

## References

- Brissett, N. C. et al. (2020). *Ku and DNA repair in bacteria*.  
- Parks, D. H. et al. (2020). Genome Taxonomy Database (GTDB).  
- UniProt Consortium (2023). UniProt reference proteomes.  
- Blum, M. et al. (2021). InterPro in 2021: improving coverage and classification.

---

## Reproducibility
All scripts, intermediate outputs, and final results are organized in this repository to ensure reproducibility of the analysis pipeline.

---

## Instructor
Dr. Saurabh Mahajan

## Author
Suraj Rathi

---
