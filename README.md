# nhej-repair-genome-analysis
Large-scale analysis of NHEJ repair systems and their association with GC content across bacterial genomes.
# 🧬 GC Content and Evolution of NHEJ Repair Systems in Bacteria

## 📌 Overview
This project investigates the distribution and evolution of bacterial Non-Homologous End Joining (NHEJ) DNA repair systems and their association with genomic GC content. The study integrates large-scale genomic datasets with multiple computational approaches to identify repair domains, map them onto phylogeny, and analyze compositional trends.

---

## 🎯 Research Objective
The primary goal of this study is to understand the biological and evolutionary factors contributing to GC content diversity in bacterial genomes, with a particular focus on the role of DNA repair systems such as NHEJ.

---

## 🔬 Key Questions
- Is the presence of NHEJ components (Ku and LigD) associated with GC-rich genomes?
- How are DNA repair systems distributed across bacterial phylogeny?
- Can large-scale datasets provide more robust insights compared to previous studies?

---

## 📊 Methodology Overview

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

## 📁 Repository Structure
nhej-repair-genome-analysis/
│
├── data/ # Raw and processed datasets
│
├── scripts/ # Analysis scripts
│ ├── hmm/
│ ├── interpro/
│ ├── phylogeny/
│ └── gc_analysis/
│
├── results/ # Output files and domain results
│ ├── hmm_profiles/
│ ├── hmmer/
│ ├── interpro/
│ ├── phylogeny/
│ └── gc/
│
├── figures/ # Final plots and visualizations
│
├── notebooks/ # Jupyter notebooks (if any)
│
├── docs/ # Supporting documentation
│
├── README.md

## 📚 References

- Brissett, N. C. et al. (2020). *Ku and DNA repair in bacteria*.  
- Parks, D. H. et al. (2020). Genome Taxonomy Database (GTDB).  
- UniProt Consortium (2023). UniProt reference proteomes.  
- Blum, M. et al. (2021). InterPro in 2021: improving coverage and classification.

## 🚀 Reproducibility
All scripts, intermediate outputs, and final results are organized in this repository to ensure reproducibility of the analysis pipeline.

---
## Instructor
Dr. Saurabh Mahajan

## Author
Suraj Rathi


---
