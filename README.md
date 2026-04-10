# Tabular ML for Somatic Variant Classification

## Overview

This repository implements a machine learning pipeline for somatic variant classification using tabular features derived from variant call data sets. The  objective is to distinguish somatic vs germline variants in tumour-only sequencing data and to derive downstream metrics such as tumour mutational burden (TMB).

The pipeline operates on annotated VCF files, generated using the mulo-wesml pipeline (https://github.com/lt11/mulo-wesml), ensuring consistency between upstream variant calling and downstream machine learning classification.

Training is performed using matched tumour–normal variant calls, enabling supervised learning in a setting where ground truth is otherwise unavailable.

The pipeline performs:

- VCF parsing and annotation extraction
- Feature engineering from variant-level, sequence-context, and functional annotations
- One-hot encoding of categorical features
- Biologically motivated filtering of germline variants
- Training and evaluation of multiple ML models
- Prediction across validation and test cohorts

---

## Input Data

The input VCF data can be prepared using the script "aux/prep-data.sh" after setting the folder paths in the "settings" section. 

A table with metadata (e.g. cancer type, sample identifier, sample ancestry), e.g. "tab/metadata.csv", is also needed.

---

## Detailed Description of the Pipeline

### Features

Key features extracted from VCF files include:

- Quantitative: DP, AF, AD, CNT
- Sequence context: substitution type, trinucleotide context
- Functional annotation: ontology categories (missense, nonsense, indels, non-coding)
- Population frequency: pop_max from dbNSFP
- Target: somatic (1) vs germline (0)

### Models

The repository supports multiple machine learning approaches:

- Gradient boosting (XGBoost, LightGBM)
- Deep tabular learning (TabNet)

These models are designed to capture nonlinear relationships and interactions between genomic features.

### Preprocessing

Categorical variables, including substitution type, trinucleotide context, and functional ontology, are transformed using one-hot encoding to produce a machine-learning–compatible feature space. The encoding scheme is fitted on the training dataset and consistently reused across validation and test data sets to ensure feature alignment. In addition, biologically motivated filtering steps are applied to reduce noise: common germline variants are removed based on population frequency (pop_max ≥ 0.01), and the remaining germline variants are further restricted to a subset of functionally relevant annotation classes.

---

## Outputs

The pipeline generates:

- Processed datasets (*-df.csv, *-df-flt.csv)
- Model files
- Prediction tables (preds-*.csv)
- ROC curves and AUC scores
- Feature importance plots

Additional tables and plots can be generated running "aux/run-r.sh". These include:
- Tumour Mutational Burden (TMB) tables
- TMB box plots
- TMB regression against the literature
- scores (precision, recall, ...) for each variant type
- scores (precision, recall, ...) for each test data set
- scores (precision, recall, ...) for each sample

---

## Installation

Clone the repository:

```
git clone https://github.com/lt11/tab-unmatched.git
cd tab-unmatched
```

Then, install the following Python modules:

- numpy 1.16.4
- pandas
- matplotlib
- seaborn
- scikit-learn
- xgboost
- lightgbm
- torch
- torchvision
- category_encoders
- joblib

---

## Usage

Once you have prepared the VCF data with "prep-data.sh", running the script is as simple as:

```
python unc-normal-drop-5.py
```

Users can flexibly tailor the feature space by editing the "optional_drop_features" parameter defined at the beginning of the script. By specifying a list of feature names (e.g. ['CNT', 'DP', 'ontology_missense']), these variables will be excluded from the training and prediction matrices prior to model fitting. If set to ['none'], no additional features are removed beyond the default metadata columns. This design enables straightforward experimentation with feature selection strategies, allowing users to assess the impact of specific variables on model performance without modifying the core pipeline logic.

---

## Citation

If you use this code, please cite:

Tattini, L., Yan, Y., Chaturvedi, N., & Appuswamy, R. (2025). Accurate Variant Classification in Tumour-Only Genomic Data Using Interpretable Tabular Models. bioRxiv, 2025-12.
doi: https://doi.org/10.64898/2025.12.09.693348
