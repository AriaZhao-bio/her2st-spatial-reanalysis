# HER2ST Breast Cancer Spatial Transcriptomics Reanalysis

## Overview

This project reanalyzes the public HER2ST breast cancer spatial transcriptomics dataset. The analysis covers 36 tissue sections from 8 patients and retains 13,302 spots after quality control. It examines ERBB2 expression, section-associated variation, and four spatial gene programs. Pathology-region comparisons use one annotated section per patient (8 sections in total).

## Data sources

- The expression matrices, H&E images, spot coordinates, and pathology annotations come from the public [HER2ST dataset on Zenodo](https://doi.org/10.5281/zenodo.4751624). The original data files are downloaded locally and are not stored in this Git repository.
- The E2F, hypoxia, and epithelial–mesenchymal transition Hallmark gene sets come from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb). Their use and redistribution are subject to the [MSigDB license terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).
- The nine-gene HER2/17q12 list and the gene-set file sources are recorded in `results/tables/gene_program_gene_set_manifest.csv`.

## Analysis workflow

1. Inspect the 36 sections, match expression spots to spatial coordinates, and perform quality control. Retain 13,302 of 13,620 spots.
2. Build one Seurat object per patient. Log-normalize expression, select variable genes, scale the data, and run PCA within each patient.
3. Examine section-associated variation and map normalized ERBB2 expression.
4. Match pathology labels to 8 annotated sections, then compare invasive-cancer and nonmalignant regions using patient-level summaries.
5. Calculate four gene-program scores and test patient-level regional differences, applying BH correction across the four programs.

## Key findings

- Section-associated variation across the first 10 PCs differed between patients, ranging from approximately 0.013 to 0.273.
- In the 8 pathology-annotated sections, the patient-level median normalized ERBB2 expression was higher in invasive-cancer regions than in pooled nonmalignant regions for all 8 patients. The median within-patient difference was 0.456 (two-sided exact Wilcoxon signed-rank test, p = 0.0078125).
- HER2/17q12, E2F proliferation, and hypoxia gene-program scores were higher in invasive-cancer regions in all 8 patients (BH-adjusted FDR = 0.0104 for each program).
- The EMT/stromal score was higher in 5 of 8 patients (BH-adjusted FDR = 0.25), providing no consistent evidence of an increase across patients.

## Limitations

- Pathology-region comparisons use only one annotated section from each of 8 patients; the other 28 sections do not have these region labels.
- Spots within the same section are correlated. Statistical tests therefore use patient-level summaries (n = 8), not individual spots as independent replicates.
- Gene-program scores are based on expression standardized within each patient. They describe relative expression patterns and do not directly measure DNA amplification, cell proliferation, or tissue oxygen levels.
- The HER2/17q12 program includes ERBB2, so its agreement with the single-gene result is not an independent validation. The EMT/stromal score may also reflect differences in tissue composition.


## Reproducing the analysis

- Tested with R 4.5.2 and Seurat 5.5.1. The scripts also use Matrix, ggplot2, patchwork, and jpeg.
- Open `her2st-spatial-reanalysis.Rproj` in RStudio and run the scripts in `scripts/` in numerical order, from `00_setup_project.R` through `14_gene_program_scores.R`.
- `01_download_data.R` and `05_download_extract_images.R` download the public files from Zenodo. Before running the extraction scripts on another computer, update the `seven_zip` path in `02_extract_data.R` and `05_download_extract_images.R` to point to a local 7-Zip executable.
- Downloaded and processed data are excluded from Git. The four small gene-set files needed for scoring are included under `data/raw/gene_sets/`.
- Summary tables and figures are written to `results/tables/` and `results/figures/`.