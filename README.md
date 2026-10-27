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

The input VCF files are those produced by the mulo-wesml pipeline, one per tumour sample, named after its fastq id. Two tables describe them:

- "tab/metadata.csv": the metadata of each sample (e.g. cancer type, sample identifier, sample ancestry); its "split" column assigns every tumour sample to one set: train, validation, test-mixtcga, test-mixanc, test-metmel or test-ffpe-dlbcl
- "tab/index-paths.txt": the path of every VCF, one per line, relative to the mulo-wesml folder, e.g. "run-9/var-calls/gt/9fb63209f7f145db94e3f57f4e510a4f.vcf.gz"

The superrunner (see Usage) checks both tables and prepares the VCFs of each set.

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

To follow the approach of the original paper, these filters are applied to the training, validation, and test sets using the ground-truth somatic/germline labels derived from matched tumour–normal calls. Both filters remove only variants labelled germline (`target == 0`); somatic variants are retained regardless of these filtering criteria.

Consequently, validation and test metrics describe performance on a population selected using ground-truth labels. They should not be interpreted as performance on an unfiltered tumour-only VCF, because the same filtering cannot be reproduced for unlabelled samples whose somatic/germline status is unknown. This label-dependent evaluation is retained for comparability with the original paper and should be considered when interpreting the reported results.

---

## Outputs

The pipeline generates:

- Processed datasets (*-df.csv, *-df-flt.csv)
- Model files
- Prediction tables (preds-*.csv)
- ROC curves and AUC scores (except for train data sets of the three ML models)
- Feature importance plots
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
- pytorch-tabnet
- category_encoders
- joblib

and the following R packages:

- data.table
- ggplot2
- this.path
- scriptName

---

## Usage

### Feature progression with the superrunner

The file "scr/superrunner.sh" runs a whole feature progression in one go. Run 1 keeps every feature, run 2 drops the first feature of the "drop_order" list in the settings of the superrunner, run 3 the first two, and so on. Every run has five replicates (a to e), each with its own seed, and the same seeds are used for every run. Finally, the superrunner collects the stats with "scr/get-stats.sh" and plots them with "scr/plot-featprog-reps.R".

Each experiment lives in its own folder, called $dir_base below:

```
$dir_base/
├── scr/                   all the code: scr/* of this repository
├── dat/
│   ├── tab/               metadata.csv and index-paths.txt
│   ├── splits/            the sample lists of each set and their checks
│   ├── vcf/               the VCFs of each set, copied and unpacked once
│   └── feat/              the feature tables, built once and shared by every run
├── run-1a/ ... run-1e/    every feature
├── run-2a/ ... run-2e/    every feature but the first of drop_order, and so on
│   └── res/ plots/ log/   the results, the TMB plots and the logs of the run
└── res/
    ├── accu-prog/         stats of the test splits, and of train and validation
    └── accu-prog-plots/   one plot per metric
```

To set up and launch an experiment:

```
mkdir -p $dir_base/scr $dir_base/dat/tab
cp scr/*.py scr/*.R scr/*.sh $dir_base/scr/
cp tab/metadata.csv tab/index-paths.txt $dir_base/dat/tab/
cd $dir_base
nohup bash scr/superrunner.sh > superrunner.log 2>&1 &
```

The settings at the top of the superrunner hold the location of the mulo-wesml folder ("vcf_src_candidates", the first that exists is used), the replicates and their seeds, and "drop_order".

The superrunner stops at the first failure, with a message pointing to the log of the step that failed. Launched again, it keeps what is already done (the data and the finished runs) and redoes an unfinished run from scratch. It never mixes results: a finished run made with other settings, or data prepared from another metadata.csv or index-paths.txt, stops it with a message saying what to remove.

### A single run

The Python script can also be run on its own. With no paths given, it reads train/, validation/, test/<subtype>/ and tab/metadata.csv in the parent of the working directory, and writes its results in res/ there:

```
python unc-normal-drop-10.py <seed> [<run_logit>] [<run_hf>] [--drop-features <feature> ...]
```

The options "--dir-vcf", "--file-meta", "--dir-out" and "--dir-feat" set the paths explicitly. With "--stage prep" the script only builds the feature tables, and with "--stage train" it only trains the models on tables built before.

Users can tailor the feature space with the "--drop-features" option. The features listed (e.g. CNT DP ontology_missense, any of "valid_drop_features" at the beginning of the script) are excluded from the training and prediction matrices prior to model fitting. With "none", the default, no features are removed beyond the default metadata columns, and an unknown feature name stops the script. This design enables straightforward experimentation with feature selection strategies, allowing users to assess the impact of specific variables on model performance without modifying the core pipeline logic.

---

## Citation

If you use this code, please cite:

Tattini, L., Yan, Y., Chaturvedi, N., & Appuswamy, R. (2025). Accurate Variant Classification in Tumour-Only Genomic Data Using Interpretable Tabular Models. bioRxiv, 2025-12.
doi: https://doi.org/10.64898/2025.12.09.693348
