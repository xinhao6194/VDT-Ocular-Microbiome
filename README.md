# VDT-Ocular-Microbiome
Custom bioinformatics and statistical scripts for the metagenomic analysis of ocular surface microbiome under prolonged VDT exposure.
# VDT-Ocular-Microbiome

Custom bioinformatics and statistical scripts for the metagenomic analysis of the ocular surface microbiome under prolonged Visual Display Terminal (VDT) exposure.

## 📖 Project Description
This repository contains the code and processed data supporting the manuscript: **"Impact of Prolonged Visual Display Terminal Usage on Ocular Surface Microbiome: A Metagenomic Analysis"**.

In accordance with Open Science and FAIR data principles, we have provided all scripts and de-identified data matrices required to reproduce the statistical analyses and figures presented in the study.

## 📁 Repository Structure

- **`scripts/`**: R scripts for downstream data analysis and visualization.
  - `alpha_diversity.R`: Script for calculating and plotting alpha diversity indices.
  - `permanova.R`: Script for PCoA and PERMANOVA beta-diversity analysis (includes `set.seed(12345)`).
  - `phylum.R`: Script for generating phylum-level relative abundance plots.
  - `genus.R`: Script for generating genus-level relative abundance plots and LEfSe-related visualizations.
  - `maaslin2_clinical_association.R`: Script for single-predictor linear mixed-effects models associating taxa with clinical metrics.

- **`data/`**: Processed data matrices and metadata.
  - `clinical_metadata_deidentified.txt`: De-identified clinical parameters (OSDI, TBUT, CFS, Meibomian Gland Score) for all 48 samples.
  - `group.txt`: Sample grouping metadata (0W, 1W, 3W).
  - `merged_abundance_table.txt`: Merged taxonomic relative abundance matrix.
  - `keggpath_abundance.community.txt`: KEGG pathway-level abundance matrix.
  - `keggmodule_abundance.community.txt`: KEGG module-level abundance matrix.

*(Note: The upstream sequencing data processing pipeline (QC, assembly, gene prediction) is referenced in the manuscript. Raw sequencing data is available at NCBI SRA under BioProject PRJNA1444215.)*

## 💻 System Requirements & Reproducibility

### Statistical Analysis (R)
The scripts were executed using R (v4.x.x). The exact versions of all packages (e.g., `vegan`, `MaAsLin2`, `ggplot2`, `dplyr`) and dependency trees are documented in `sessionInfo_summary.txt`.

### Key Statistical Settings
- **Random Seed:** To ensure exact reproducibility of permutation-based analyses, `set.seed(12345)` was explicitly set prior to PERMANOVA (`adonis2` in the `vegan` package).
- **MaAsLin2 Models:** To avoid collinearity among clinical variables, each metric was fitted in a separate single-predictor linear mixed-effects model: `Feature ~ Clinical_Metric + (1|Subject)`. Multiple testing correction used the Benjamini-Hochberg procedure.

## 🚀 How to Run
1. Open the R scripts in the `scripts/` folder.
2. Ensure the working directory is set to the root of this repository.
3. Run the scripts sequentially to reproduce the plots and statistical outputs.

## ✉️ Contact
For any questions regarding the code or data, please contact:
- **Ling Li**: Liling498@126.com
- **Wei Chen**: chenweimd@wmu.edu.cn
