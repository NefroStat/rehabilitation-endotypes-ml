# Machine learning reveals systemic endotypes of extracellular matrix remodeling and endothelial dysfunction independent of injury severity in delayed combat trauma

## Overview
This repository contains the computational pipeline and analytical code for the study evaluating subclinical pathophysiological endotypes in combat trauma survivors during the delayed rehabilitation phase. The methodology integrates classical non-parametric statistics, dimensionality reduction, and non-linear machine learning algorithms to assess the interaction between extracellular matrix (ECM) turnover and endothelial dysfunction. 

## Data Availability
In compliance with the ethical approval granted by the Bioethics Committee of Danylo Halytsky Lviv National Medical University (Protocol No. 9) and regulations concerning military personnel data confidentiality, the raw clinical dataset is not publicly available in this repository. 
* To ensure methodological transparency, a detailed specification of the required variables is provided in `data/DATA_DICTIONARY.md`.
* Researchers seeking to reproduce the analysis must request access to the de-identified dataset via the corresponding author, subject to a formal data-sharing agreement.

## System Requirements
The analysis was performed in the R programming environment.
* **R version:** 4.5.3 or higher.
* **Core Dependencies:** `dplyr`, `gtsummary`, `cluster`, `pwr`, `robustbase`, `glmnet`, `randomForest`, `xgboost`, `caret`, `psych`, `factoextra`, `ggplot2`, `patchwork`, `openxlsx`, `tidyr`, `scales`.

Dependencies can be installed using the standard `install.packages()` function. Ensure the `renv.lock` file is utilized if a strict environment replication is required.

## Repository Architecture

├── README.md
├── LICENSE
├── .gitignore
├── renv.lock
├── data/
│   └── DATA_DICTIONARY.md           # Metadata and variable specifications
├── scripts/
│   ├── 01_data_preprocessing.R      # Data cleaning and index calculation
│   ├── 02_statistical_analysis.R    # NISS stratification and non-parametric tests
│   ├── 03_hierarchical_clustering.R # Ward's clustering and endotype assignment
│   ├── 04_pca_dimensionality_reduction.R # PCA, Robust Regression, and LASSO
│   ├── 05_xgboost_shap_modeling.R   # Random Forest, XGBoost, and SHAP extraction
│   └── 06_data_visualization.R      # Publication-ready figure generation
└── output/
    ├── figures/                     # Exported TIFF plots (300 DPI)
    └── tables/                      # Comprehensive Excel statistical outputs


## Execution Protocol
To execute the pipeline, a dataset formatted exactly according to `DATA_DICTIONARY.md` must be placed in the `data/` directory and named `markers_results_comma.csv`. 

The scripts must be executed sequentially to maintain data dependency logic:
1. Run `01_data_preprocessing.R` to format variables and export `preprocessed_data.rds`.
2. Run `02_statistical_analysis.R` to compute anatomical severity differences and export `niss_stats.rds`.
3. Run `03_hierarchical_clustering.R` to algorithmically define Endotype 1 and Endotype 2, exporting `clustered_data.rds`.
4. Run `04_pca_dimensionality_reduction.R` to compute eigenvectors and linear matrices, exporting `pca_linear_data.rds`.
5. Run `05_xgboost_shap_modeling.R` to evaluate non-linear feature importance. This script generates the final statistical workbook `Final_Analysis_Complete.xlsx` in the `output/tables/` directory.
6. Run `06_data_visualization.R` to render and save Figures 1-4 to the `output/figures/` directory.

## License
The code in this repository is licensed under the MIT License. See the `LICENSE` file for details.

## Citation
If you utilize this code or methodology in your research, please cite the original article and the Zenodo repository DOI.
* [Placeholder for Article Citation]
* [Placeholder for Zenodo DOI]