from __future__ import division

## settings ------------------------------------------------------------------

### 'none' or any valid feature e.g. 'CNT', 'ontology_1', 'DP'
### (of type list)
optional_drop_features = ['none']
### important features
# optional_drop_features = ['CNT']
# optional_drop_features = ['CNT',
#                           'ontology_NR',
#                           'DP',
#                           'ontology_inframe_indel', 
#                           'AF',
#                           'ontology_frame_shift_indel',
#                           'AD',
#                           'pop_max',
#                           'subs_type_NA',
#                           'ontology_nonsense',
#                           'subs_type_G>A',
#                           'trinucleotide_TCA',
#                           'trinucleotide_GGG',
#                           'subs_type_G>T',
#                           'subs_type_C>T',
#                           'ontology_missense',
#                           'subs_type_G>C',
#                           'trinucleotide_TCC',
#                           'trinucleotide_GGA',
#                           'trinucleotide_TGA',
#                           'trinucleotide_GTG',  
#                           'subs_type_C>A']
### useless features
# optional_drop_features = ['trinucleotide_GAC']
# optional_drop_features = ['trinucleotide_GAC', 'trinucleotide_ACA',
#                           'trinucleotide_GCC', 'subs_type_A>T']

## common methods -------------------------------------------------------------

def get_file_names(folder_path):
    file_names = []
    for root, dirs, files in os.walk(folder_path):
        for file_name in files:
            if "tbi" not in file_name and "vcf" in file_name:
                file_names.append(os.path.join(root, file_name))
    return file_names

### accuracy method
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score as ras

def accuracy(model, test_val, preds):
    print("Area under the precision-recall curve ", ras(test_val.values, preds))
    fpr, tpr, _ = roc_curve(test_val.values, preds)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, 'blue', label = 'ROC curve (area = %0.2f)' % roc_auc)
    plt.scatter(fpr, tpr, color='red', s=10)
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(model + ' ROC curve')
    plt.legend(loc="lower right")
    plt.savefig(output_dir + "/roc-" + model + ".pdf")
    plt.show()

### plot feature importance
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
### all features
def myplot_importance(feature_names, feature_importance, model_name,
                      output_dir):
    ### make a DataFrame for better visualization
    importance_df = pd.DataFrame({'Feature': feature_names, 'Importance':
                                  feature_importance})
    ### sort the DataFrame by importance
    importance_df = importance_df.sort_values(by='Importance', ascending=False)
    ### plot the feature importance
    plt.figure(figsize=(10, 8))
    sns.barplot(x='Importance', y='Feature', data=importance_df, color='black')
    plt.xlabel('Importance Score (a.u.)', fontsize=22)
    plt.ylabel('Features', fontsize=22)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(output_dir + '/importance-' + model_name + '.pdf')
    plt.show()

### number of top features to plot
n_top_feat = 20
def myplot_importance_top(feature_names, feature_importance, model_name,
                          output_dir):
    ### make a DataFrame for better visualization
    importance_df = pd.DataFrame({'Feature': feature_names, 'Importance':
                                  feature_importance})
    ### sort the DataFrame by importance
    importance_df = importance_df.sort_values(by='Importance', ascending=False)
    ### plot the feature importance
    plt.figure(figsize=(10, 8))
    sns.barplot(x='Importance', y='Feature', data=importance_df, color='black')
    plt.xlabel('Importance Score (a.u.)', fontsize=22)
    plt.ylabel('Features', fontsize=22)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=22)    
    plt.tight_layout()
    plt.savefig(output_dir + '/importance-' + model_name + '.pdf')
    plt.show()

### load the vcf files
def load_variants(input_dir):
    file_names = get_file_names(input_dir)
    print('Nb of VCF files:', len(file_names))
    dfs = []
    for file_name in file_names:
        with open(file_name, 'r') as file:
            for line_number, line in enumerate(file):
                if line.startswith("#CHROM"):
                    break
        df = pd.read_csv(file_name, sep="\t", skiprows=line_number,
                         quoting=csv.QUOTE_NONE)
        df.rename(columns={df.columns[9]: 'SAMPLE'}, inplace=True)
        my_name = os.path.basename(file_name)
        fastq_id = re.sub(r'\.vcf$', '', my_name)
        df['fastq id'] = fastq_id
        dfs.append(df)
    raw = pd.concat(dfs, ignore_index=True)
    raw = pd.merge(raw, info[['fastq id', 'subtype']], on='fastq id',
                   how='left')
    print('Nb of total variants: ', len(raw))
    return raw

import re

### extraction of the string after key until the first occurrence of ';'
def extract_key(text, key):
    pattern = key + '=(.*?);'
    result = text.str.extract(pattern, expand=False)
    ### replace NaNs with empty string
    return result.fillna('')

def fill_target(row):
    if row['STATUS'] == 'Somatic':
        return 1
    elif row['STATUS'] == 'Germline':
        return 0

def fill_cnt(row):
    if row['CNT'] == '.':
        return 0
    else:
        return int(row['CNT'])

def fill_substitution(row):
    if row['TYPE'] == 'SNV':
        return row['REF'] + '>' + row['ALT']
    else:
        return 'NA'

def fill_trinucleotide(row):
    if row['TYPE'] == 'SNV':
        return row['LSEQ'][-1] + row['REF'] + row['RSEQ'][0]
    else:
        return 'NA'

### feature ontology
ontology_type = {
    'missense':[
        'stop_lost',
        'start_lost',
        'missense_variant',
        'missense_variant&splice_region_variant'],
    'nonsense':[
        'stop_gained',
        'frameshift_variant&stop_gained',
        'stop_gained&splice_region_variant'],
    'inframe_indel':[
        'disruptive_inframe_deletion',
        'conservative_inframe_deletion',
        'conservative_inframe_insertion'],
    'frame_shift_indel':[
        'frameshift_variant&splice_donor_variant&' ### yep, no comma
        'splice_region_variant&intron_variant',
        'frameshift_variant&splice_region_variant',
        'frameshift_variant&stop_lost&splice_region_variant',
        'frameshift_variant'],
    'NR':[
        'upstream_gene_variant',
        'intron_variant',
        'splice_region_variant&intron_variant',
        '3_prime_UTR_variant',
        'synonymous_variant',
        '5_prime_UTR_variant',
        'downstream_gene_variant',
        'splice_region_variant&synonymous_variant',
        'splice_acceptor_variant&intron_variant',
        'splice_region_variant&non_coding_transcript_exon_variant',
        'non_coding_transcript_exon_variant',
        'intergenic_region',
        '5_prime_UTR_premature_start_codon_gain_variant',
        'initiator_codon_variant',
        'splice_donor_variant&intron_variant',
        'splice_region_variant',
        'splice_region_variant&stop_retained_variant',
        'stop_retained_variant',
        'splice_donor_variant&conservative_inframe_deletion&' ### yep, no comma
        'splice_region_variant&intron_variant',
        'splice_acceptor_variant&conservative_inframe_deletion&'
        'splice_region_variant&intron_variant']}

def check_ontology(val):
    for sublist_name, sublist in ontology_type.items():
        if val in sublist:
            return sublist_name
    return 'NA'

### feature extraction for pop_max
dbNSFP = ['dbNSFP_ExAC_NFE_AF',
          'dbNSFP_ExAC_SAS_AF',
          'dbNSFP_ExAC_Adj_AF',
          'dbNSFP_1000Gp3_AMR_AF',
          'dbNSFP_ExAC_AFR_AF',
          'dbNSFP_1000Gp3_AF',
          'dbNSFP_1000Gp3_EAS_AF',
          'dbNSFP_ExAC_AF',
          'dbNSFP_ExAC_FIN_AF',
          'dbNSFP_1000Gp3_EUR_AF',
          'dbNSFP_ExAC_AMR_AF',
          'dbNSFP_1000Gp3_AFR_AF',
          'dbNSFP_ESP6500_AA_AF',
          'dbNSFP_1000Gp3_SAS_AF',
          'dbNSFP_ExAC_EAS_AF',
          'dbNSFP_ESP6500_EA_AF']

def find_dbNSFP(string):
    pattern = r'dbNSFP_[^;]+;'
    matches = re.findall(pattern, string)
    return matches

def get_pop_max(dbs):
    val = []
    for db in dbs:
        db_name = db.split('=')[0]
        if db_name in dbNSFP:
            db_value = db.split('=')[1]
            if db_value[-1] == ';':
                db_value = db_value[:-1]
                if ',' in db_value:
                    db_value = db_value.split(',')[0]
                val.append(db_value)
    if len(val) == 0:
        return 0
    else:
        float_values = [float(x) for x in val]
        return max(float_values)

from category_encoders import OneHotEncoder
import joblib

AD_index = 3
BIAS_index = 7

def extract_feature(raw, mode, output_dir):
    raw['TYPE'] = extract_key(raw['INFO'], 'TYPE')
    raw['DP'] = extract_key(raw['INFO'], 'DP').astype(int)
    raw['AF'] = extract_key(raw['INFO'], 'AF').astype(float)
    raw['AD'] = raw['SAMPLE'].str.split(':').str[AD_index]
    raw['AD'] = raw['AD'].apply(lambda x: list(map(int, x.split(','))))
    raw['AD_max'] = raw['AD'].apply(max)
    raw['AD_sum'] = raw['AD'].apply(sum)
    raw['AD'] = raw['AD_max']/raw['AD_sum']
    raw['LSEQ'] = extract_key(raw['INFO'], 'LSEQ')
    raw['RSEQ'] = extract_key(raw['INFO'], 'RSEQ')
    raw['STATUS'] = extract_key(raw['INFO'], 'STATUS')
    raw['SAMPLE_ID'] = extract_key(raw['INFO'], 'SAMPLE')
    raw['nb_variants'] = raw.groupby('SAMPLE_ID')['SAMPLE_ID'].transform('count')
    ### convert somatic and germline to 1 and 0, respectively
    raw['STATUS'] = raw['STATUS'].apply(lambda x: 'Somatic' if x in
                                        ['StrongSomatic',
                                         'LikelySomatic'] else x)
    raw['STATUS'] = raw['STATUS'].apply(lambda x: 'Germline' if x == '' else x)
    raw['target'] = raw.apply(fill_target, axis=1)
    ### get CNT
    raw['CNT'] = extract_key(raw['INFO'], 'CNT').str.split(',').str[0]
    raw['CNT'] = raw.apply(fill_cnt, axis=1)
    raw['subs_type'] = raw.apply(fill_substitution, axis=1)
    raw['trinucleotide'] = raw.apply(fill_trinucleotide, axis=1)
    raw['ANN'] = extract_key(raw['INFO'], 'ANN').str.split('|').str[1]
    raw['ontology'] = raw['ANN'].apply(check_ontology)
    ### get pop_max
    raw['raw_pop_max'] =  raw['INFO'].apply(find_dbNSFP)
    raw['pop_max'] = raw['raw_pop_max'].apply(get_pop_max)
    ### TODO
    ### raw['BIAS'] = raw['SAMPLE'].str.split(':').str[BIAS_index]
    ### raw['BIAS_REF'] = raw['BIAS'].str.split(',').str[0].astype(int)
    ### raw['BIAS_ALT'] = raw['BIAS'].str.split(',').str[1].astype(int)
    ### feature study 
    print('--- Type of variants quality [PASS or NOT]---')
    print(raw['FILTER'].unique())
    print('--- STATUS of varints ---')
    print(raw['target'].unique())
    print('--- TYPE of varints ---')
    print(raw['TYPE'].unique())
    print('--- Number of varints in each sample---')
    print(raw['nb_variants'].unique())
    print('--- ontology possible values [check there is no '']---')
    print(raw['ontology'].unique())
    ### print('--- BIAS all possible values---')
    ### print(raw['BIAS'].unique())
    ## distribution of varints type germline/somatic ---------------------------
    value_counts = raw['STATUS'].value_counts()
    plt.figure(figsize=(4, 4))
    plt.pie(value_counts, labels=value_counts.index,
            autopct='%1.1f%%', startangle=140)
    plt.title('Pie Chart of variants type')
    ### equal aspect ratio ensures that pie is drawn as a circle
    plt.axis('equal')
    plt.savefig(output_dir + '/' + mode + '-variants-pie-chart.png')
    plt.show()
    ### output formated data
    col = ['DP', 'AF', 'AD', 'CNT', 'subs_type', 'trinucleotide',
           'nb_variants', 'ontology', 'pop_max', 'target', 'subtype',
           'fastq id', 'ANN', '#CHROM', 'POS']
    data = raw[col].copy()
    print(data.dtypes)
    data.to_csv(output_dir + '/' + mode + '-formatted-variants.csv',
                index=False)
    ### output OneHotEncoded data
    encoder_path = output_dir + '/onehot_encoder.pkl'
    if mode == 'train':
        ### replaces categorical columns e.g. subs_type 
        ### with multiple one-hot encoded columns e.g. subs_type_A>T
        encoder = OneHotEncoder(cols=['subs_type', 'trinucleotide', 'ontology'],
                                use_cat_names=True)
        data = encoder.fit_transform(data)
        joblib.dump(encoder, encoder_path)
    else:
        encoder = joblib.load(encoder_path)
        data = encoder.transform(data)
    train_X = data.drop(columns='target')
    train_Y = data['target']
    print("set shape: ")
    print(data.shape)
    print("somatic variants Nb: " + str(len(train_Y[train_Y==1])))
    print("germline variants Nb: " + str(len(train_Y[train_Y==0])))
    data.to_csv(output_dir + '/' + mode + '-df.csv', index=False)
    return data

## common data -----------------------------------------------------------------

import pandas as pd
import sys, csv, os, re
import matplotlib.pyplot as plt
import numpy as np
import os

current_directory = os.getcwd()
workspace = os.path.dirname(current_directory) + '/'

output_dir = workspace + 'res/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

### load info
info_dir = workspace + 'tab/metadata.csv'
info = pd.read_csv(info_dir, sep=",", quotechar='"', engine='python')
info.head(2)

### filtering conditions
ann_valid = {"stop_lost",
             "start_lost",
             "missense_variant",
             "missense_variant&splice_region_variant",
             "stop_gained",
             "frameshift_variant&stop_gained",
             "stop_gained&splice_region_variant",
             "disruptive_inframe_deletion",
             "conservative_inframe_deletion",
             "conservative_inframe_insertion",
             "frameshift_variant&splice_donor_variant" ### yep, no comma
             "&splice_region_variant&intron_variant",
             "frameshift_variant&splice_region_variant",
             "frameshift_variant&stop_lost&splice_region_variant",
             "frameshift_variant"}

## training data preparation ---------------------------------------------------

mode = 'train'
print("processing: " + mode)
input_dir = workspace + mode + '/'

### load vcf
raw = load_variants(input_dir)
data = extract_feature(raw, mode, output_dir)
data.drop('nb_variants', axis=1, inplace=True)
### filter common germline variants
data.drop(data[(data['pop_max'] >= 0.01) &
               (data['target'] == 0)].index,
               inplace=True)
count = (data['target'] == 0).sum()
print("post pop_max filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post pop_max filters somatic variants Nb: " + str(count))
ind_out = data[(data['target'] == 0) & (~data['ANN'].isin(ann_valid))].index
data.drop(index=ind_out, inplace=True)
count = (data['target'] == 0).sum()
print("post ANN filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post ANN filters somatic variants Nb: " + str(count))
data.to_csv(output_dir + '/' + mode + '-df-flt.csv', index = False)
data.head(2)

## validation data preparation -------------------------------------------------

mode = 'validation'
print("processing: " + mode)
input_dir = workspace + mode + '/'
output_dir = workspace + 'res/'

### load vcf
raw = load_variants(input_dir)
data = extract_feature(raw, mode, output_dir)
data.drop('nb_variants', axis=1, inplace=True)
### filter common germline variants
data.drop(data[(data['pop_max'] >= 0.01) &
               (data['target'] == 0)].index,
               inplace=True)
count = (data['target'] == 0).sum()
print("post pop_max filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post pop_max filters somatic variants Nb: " + str(count))
ind_out = data[(data['target'] == 0) & (~data['ANN'].isin(ann_valid))].index
data.drop(index=ind_out, inplace=True)
count = (data['target'] == 0).sum()
print("post ANN filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post ANN filters somatic variants Nb: " + str(count))
data.to_csv(output_dir + '/' + mode + '-df-flt.csv', index = False)
data.head(2)

## test data preparation -------------------------------------------------------

subtype = 'melanoma'
mode = 'test-' + subtype
print("processing: " + mode)
input_dir = workspace + 'test/' + subtype + '/'
output_dir = workspace + 'res/'

### load vcf
raw = load_variants(input_dir)
data = extract_feature(raw, mode, output_dir)
data.drop('nb_variants', axis=1, inplace=True)
### filter common germline variants
data.drop(data[(data['pop_max'] >= 0.01) &
               (data['target'] == 0)].index,
               inplace=True)
count = (data['target'] == 0).sum()
print("post pop_max filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post pop_max filters somatic variants Nb: " + str(count))
ind_out = data[(data['target'] == 0) & (~data['ANN'].isin(ann_valid))].index
data.drop(index=ind_out, inplace=True)
count = (data['target'] == 0).sum()
print("post ANN filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post ANN filters somatic variants Nb: " + str(count))
data.to_csv(output_dir + '/' + mode + '-df-flt.csv', index = False)
data.head(2)

subtype = 'mixtcga'
mode = 'test-' + subtype
print("processing: " + mode)
input_dir = workspace + 'test/' + subtype + '/'
output_dir = workspace + 'res/'

### load vcf
raw = load_variants(input_dir)
data = extract_feature(raw, mode, output_dir)
data.drop('nb_variants', axis=1, inplace=True)
### filter common germline variants
data.drop(data[(data['pop_max'] >= 0.01) &
               (data['target'] == 0)].index,
               inplace=True)
count = (data['target'] == 0).sum()
print("post pop_max filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post pop_max filters somatic variants Nb: " + str(count))
ind_out = data[(data['target'] == 0) & (~data['ANN'].isin(ann_valid))].index
data.drop(index=ind_out, inplace=True)
count = (data['target'] == 0).sum()
print("post ANN filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post ANN filters somatic variants Nb: " + str(count))
data.to_csv(output_dir + '/' + mode + '-df-flt.csv', index = False)
data.head(2)

## predict data preparation ----------------------------------------------------

subtype = 'tnbc'
mode = 'test-' + subtype
print("processing: " + mode)
input_dir = workspace + 'test/' + subtype + '/'
output_dir = workspace + 'res/'

### load vcf
raw = load_variants(input_dir)
data = extract_feature(raw, mode, output_dir)
data.drop('nb_variants', axis=1, inplace=True)
### filter common germline variants
data.drop(data[(data['pop_max'] >= 0.01) &
               (data['target'] == 0)].index,
               inplace=True)
count = (data['target'] == 0).sum()
print("post pop_max filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post pop_max filters somatic variants Nb: " + str(count))
ind_out = data[(data['target'] == 0) & (~data['ANN'].isin(ann_valid))].index
data.drop(index=ind_out, inplace=True)
count = (data['target'] == 0).sum()
print("post ANN filters germline variants Nb: " + str(count))
count = (data['target'] == 1).sum()
print("post ANN filters somatic variants Nb: " + str(count))
data.to_csv(output_dir + '/' + mode + '-df-flt.csv', index = False)
data.head(2)

## load formatted data ---------------------------------------------------------

import pandas as pd

current_directory = os.getcwd()
workspace = os.path.dirname(current_directory) + '/'

output_dir = workspace + 'res/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

train_df = pd.read_csv(output_dir + '/train-df-flt.csv')
train_X = train_df.drop(columns='target')
train_Y = train_df['target']
train_pred = train_df.copy()

validation_df = pd.read_csv(output_dir + '/validation-df-flt.csv')
validation_X = validation_df.drop(columns='target')
validation_Y = validation_df['target']
validation_pred = validation_df.copy()

test_df_m = pd.read_csv(output_dir + '/test-melanoma-df-flt.csv')
test_X_m = test_df_m.drop(columns='target')
test_Y_m = test_df_m['target']
test_pred_m = test_df_m.copy()

test_df_mix = pd.read_csv(output_dir + '/test-mixtcga-df-flt.csv')
test_X_mix = test_df_mix.drop(columns='target')
test_Y_mix = test_df_mix['target']
test_pred_mix = test_df_mix.copy()

test_df = pd.read_csv(output_dir + '/test-tnbc-df-flt.csv')
test_X = test_df.drop(columns='target')
test_pred = test_df.copy()

base_drop_cols = ['subtype', 'fastq id', 'ANN', '#CHROM', 'POS']
if optional_drop_features == ['none']:
    drop_cols = base_drop_cols
else:
    drop_cols = base_drop_cols + optional_drop_features

### keep only columns that exist in the current DataFrame
drop_cols = [c for c in drop_cols if c in train_X.columns]
print("Dropping columns:", drop_cols)

## xgboost ---------------------------------------------------------------------

### pip install xgboost
### pip install numpy==1.16.4
import xgboost as xgb
from xgboost import XGBClassifier
from xgboost import plot_importance

### create model instance
bst = XGBClassifier(learning_rate=0.05, n_estimators=750, max_depth=30, 
                    eval_set=[(validation_X, validation_Y)],
                    objective='binary:logistic')
                    ### avoid overfit, if 10 round no improve, stop
                    ### early_stopping_rounds=10

### fit model
bst.fit(train_X.drop(columns=drop_cols), train_Y)
bst.save_model(output_dir + 'model-xgboost.json')

### get feature importance
booster = bst.get_booster()
score_dict = booster.get_score(importance_type='gain')
### build a complete dictionary with all features, filling missing ones with 0
feature_list = train_X.drop(columns=drop_cols).columns
full_score_dict = {feat: score_dict.get(feat, 0) for feat in feature_list}
df_importance = pd.DataFrame(
    list(full_score_dict.items()),
    columns=['feature', 'importance']
    ).sort_values(by='importance', ascending=False)

### plot feature importance
myplot_importance(
    feature_names=df_importance['feature'].tolist(),
    feature_importance=df_importance['importance'].tolist(),
    model_name='xgboost',
    output_dir=output_dir
)

### top feature importance
top_feat = df_importance.head(n_top_feat)
myplot_importance_top(
    feature_names=top_feat['feature'].tolist(),
    feature_importance=top_feat['importance'].tolist(),
    model_name='xgboost-top',
    output_dir=output_dir
)

# plot_importance(bst)
# plt.savefig(output_dir + '/importance-xgboost.png')
# plt.show()

### make predictions
preds_xgboost = bst.predict_proba(validation_X.drop(columns=drop_cols))[:, 1]
### accuracy check
accuracy('xgboost', validation_Y, preds_xgboost)

### train
preds_xgboost = bst.predict_proba(train_X.drop(columns=drop_cols))[:, 1]
train_pred['xgm_preds'] = np.where(preds_xgboost>0.5,'somatic','germline')
nb_germline = len(train_Y[train_Y==0])
nb_somatic = len(train_Y[train_Y==1])
nb_germline_pred = len(train_pred[train_pred['xgm_preds']=='germline'])
nb_somatic_pred = len(train_pred[train_pred['xgm_preds']=='somatic'])
print('[training set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: '
      + str(nb_somatic))
print('[training set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### validate
preds_xgboost = bst.predict_proba(validation_X.drop(columns=drop_cols))[:, 1]
validation_pred['xgm_preds'] = np.where(preds_xgboost>0.5,'somatic','germline')
nb_germline = len(validation_Y[validation_Y==0])
nb_somatic = len(validation_Y[validation_Y==1])
nb_germline_pred = len(validation_pred[validation_pred['xgm_preds']=='germline'])
nb_somatic_pred = len(validation_pred[validation_pred['xgm_preds']=='somatic'])
print('[validation set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[validation set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test melanoma
preds_xgboost = bst.predict_proba(test_X_m.drop(columns=drop_cols))[:, 1]
test_pred_m['xgm_preds'] = np.where(preds_xgboost>0.5,'somatic','germline')
nb_germline = len(test_Y_m[test_Y_m==0])
nb_somatic = len(test_Y_m[test_Y_m==1])
nb_germline_pred = len(test_pred_m[test_pred_m['xgm_preds']=='germline'])
nb_somatic_pred = len(test_pred_m[test_pred_m['xgm_preds']=='somatic'])
print('[test melanoma set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[test melanoma set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test mixtcga
preds_xgboost = bst.predict_proba(test_X_mix.drop(columns=drop_cols))[:, 1]
test_pred_mix['xgm_preds'] = np.where(preds_xgboost>0.5,'somatic','germline')
nb_germline = len(test_Y_mix[test_Y_mix==0])
nb_somatic = len(test_Y_mix[test_Y_mix==1])
nb_germline_pred = len(test_pred_mix[test_pred_mix['xgm_preds']=='germline'])
nb_somatic_pred = len(test_pred_mix[test_pred_mix['xgm_preds']=='somatic'])
print('[test mixtcga set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[test mixtcga set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test tnbc data
preds_xgboost = bst.predict_proba(test_X.drop(columns=drop_cols))[:, 1]
test_pred['xgm_preds'] = np.where(preds_xgboost>0.5,'somatic','germline')
nb_germline = len(test_pred[test_pred['xgm_preds']=='germline'])
nb_somatic = len(test_pred[test_pred['xgm_preds']=='somatic'])
print('[test TNBC set] predicted nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))

## lightgbm --------------------------------------------------------------------

### pip install lightgbm
import lightgbm as lgb
from sklearn.metrics import mean_squared_error
import seaborn as sns

params = {'objective': 'binary','metric': 'auc', 'boosting_type': 'gbdt',
          'num_leaves': 31, 'learning_rate': 0.05, 'feature_fraction': 0.9}
### fit model
gbm = lgb.train(params,
               lgb.Dataset(train_X.drop(columns=drop_cols), train_Y),
               num_boost_round=100,
               valid_sets=lgb.Dataset(validation_X.drop(columns=drop_cols),
                                      validation_Y))
gbm.save_model(output_dir+ '/model-lgbm.txt')

### plot feature importance
myplot_importance(gbm.feature_name(), gbm.feature_importance(), 'lgbm',
                  output_dir)

### create a DataFrame from lightgbm feature importances
df_importance_lgbm = pd.DataFrame({
    'feature': gbm.feature_name(),
    'importance': gbm.feature_importance()
})

### sort and select top features
top_lgbm = df_importance_lgbm.sort_values(by='importance', 
                                          ascending=False).head(n_top_feat)

### plot top features
myplot_importance_top(
    feature_names=top_lgbm['feature'].tolist(),
    feature_importance=top_lgbm['importance'].tolist(),
    model_name='lgbm-top',
    output_dir=output_dir
)

### predict
preds_lgbm = gbm.predict(validation_X.drop(columns=drop_cols),
num_iteration=gbm.best_iteration)
### accuracy
accuracy('lgbm', validation_Y, preds_lgbm)

### train
preds_lgbm = gbm.predict(train_X.drop(columns=drop_cols),
num_iteration=gbm.best_iteration)
train_pred['lgbm_preds'] = np.where(preds_lgbm>0.5,'somatic','germline')
nb_germline = len(train_Y[train_Y==0])
nb_somatic = len(train_Y[train_Y==1])
nb_germline_pred = len(train_pred[train_pred['lgbm_preds']=='germline'])
nb_somatic_pred = len(train_pred[train_pred['lgbm_preds']=='somatic'])
print('[training set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: '
      + str(nb_somatic))
print('[training set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### validate
preds_lgbm = gbm.predict(validation_X.drop(columns=drop_cols),
num_iteration=gbm.best_iteration)
validation_pred['lgbm_preds'] = np.where(preds_lgbm>0.5,'somatic','germline')
nb_germline = len(validation_Y[validation_Y==0])
nb_somatic = len(validation_Y[validation_Y==1])
nb_germline_pred = len(validation_pred[validation_pred['lgbm_preds']==
                                       'germline'])
nb_somatic_pred = len(validation_pred[validation_pred['lgbm_preds']==
                                      'somatic'])
print('[validation set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[validation set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test melanoma
preds_lgbm = gbm.predict(test_X_m.drop(columns=drop_cols),
num_iteration=gbm.best_iteration)
test_pred_m['lgbm_preds'] = np.where(preds_lgbm>0.5,'somatic','germline')
nb_germline = len(test_Y_m[test_Y_m==0])
nb_somatic = len(test_Y_m[test_Y_m==1])
nb_germline_pred = len(test_pred_m[test_pred_m['lgbm_preds']=='germline'])
nb_somatic_pred = len(test_pred_m[test_pred_m['lgbm_preds']=='somatic'])
print('[test melanoma set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[test melanoma set] predicted nb of germline: ' 
      + str(nb_germline_pred)
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test mixtcga
preds_lgbm = gbm.predict(test_X_mix.drop(columns=drop_cols),
num_iteration=gbm.best_iteration)
test_pred_mix['lgbm_preds'] = np.where(preds_lgbm>0.5,'somatic','germline')
nb_germline = len(test_Y_mix[test_Y_mix==0])
nb_somatic = len(test_Y_mix[test_Y_mix==1])
nb_germline_pred = len(test_pred_mix[test_pred_mix['lgbm_preds']=='germline'])
nb_somatic_pred = len(test_pred_mix[test_pred_mix['lgbm_preds']=='somatic'])
print('[test mixtcga set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[test mixtcga set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test tnbc data
preds_lgbm = gbm.predict(test_X.drop(columns=drop_cols),
num_iteration=gbm.best_iteration)
test_pred['lgbm_preds'] = np.where(preds_lgbm>0.5,'somatic','germline')
nb_germline = len(test_pred[test_pred['lgbm_preds']=='germline'])
nb_somatic = len(test_pred[test_pred['lgbm_preds']=='somatic'])
print('[test TNBC set] predicted nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))

## tabnet ----------------------------------------------------------------------

### pip install torch, torchvision
from pytorch_tabnet.tab_model import TabNetClassifier
import torch

X_train = train_X.drop(columns=drop_cols).to_numpy()
Y_train = train_Y.to_numpy().squeeze()
X_validation = validation_X.drop(columns=drop_cols).to_numpy()
Y_validation = validation_Y.to_numpy().squeeze()
X_test_m = test_X_m.drop(columns=drop_cols).to_numpy()
Y_test_m = test_Y_m.to_numpy().squeeze()
X_test_mix = test_X_mix.drop(columns=drop_cols).to_numpy()
Y_test_mix = test_Y_mix.to_numpy().squeeze()
X_test = test_X.drop(columns=drop_cols).to_numpy()

### classifier
classifier = TabNetClassifier(n_d=24, n_a=24, n_steps=4, gamma=1.5,
                              n_independent=2, n_shared=2,
                              lambda_sparse=1e-4, momentum=0.3,
                              clip_value=2., optimizer_fn=torch.optim.Adam,
                              scheduler_params={"gamma": 0.95, "step_size": 20},
                              scheduler_fn=torch.optim.lr_scheduler.StepLR,
                              epsilon=1e-15)

classifier.fit(X_train=X_train, y_train=Y_train,
               eval_set=[(X_train, Y_train),
                         (X_validation, Y_validation)],
                         eval_name=['train', 'valid'], 
               max_epochs=100, patience=100,
               batch_size=4000, virtual_batch_size=256,
               eval_metric=['auc'])

### model saved as tabnet_model.zip
classifier.save_model(output_dir+ '/model-tabnet') 
### feature importance
feature_importance = classifier.feature_importances_
feature_names = train_X.drop(columns=drop_cols).columns
myplot_importance(feature_names, feature_importance, 'tabnet', output_dir)

### get top features
top_tabnet = pd.DataFrame({
    'feature': feature_names,
    'importance': feature_importance
}).sort_values(by='importance', ascending=False).head(n_top_feat)

### plot top features
myplot_importance_top(
    feature_names=top_tabnet['feature'].tolist(),
    feature_importance=top_tabnet['importance'].tolist(),
    model_name='tabnet-top',
    output_dir=output_dir
)

### predict
preds_tabnet = classifier.predict_proba(X_validation)[:, 1]
### accuracy
accuracy('tabnet', validation_Y, preds_tabnet)

### train
preds = classifier.predict_proba(X_train)[:, 1] 
train_pred['tabnet_preds'] = np.where(preds>0.5,'somatic','germline')
nb_germline = len(train_Y[train_Y==0])
nb_somatic = len(train_Y[train_Y==1])
nb_germline_pred = len(train_pred[train_pred['tabnet_preds']=='germline'])
nb_somatic_pred = len(train_pred[train_pred['tabnet_preds']=='somatic'])
print('[training set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[training set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### validate
preds = classifier.predict_proba(X_validation)[:, 1] 
validation_pred['tabnet_preds'] = np.where(preds>0.5,'somatic','germline')
nb_germline = len(validation_Y[validation_Y==0])
nb_somatic = len(validation_Y[validation_Y==1])
nb_germline_pred = len(validation_pred[validation_pred['tabnet_preds']==
                                       'germline'])
nb_somatic_pred = len(validation_pred[validation_pred['tabnet_preds']==
                                      'somatic'])
print('[validation set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[validation set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test melanoma
preds = classifier.predict_proba(X_test_m)[:, 1] 
test_pred_m['tabnet_preds'] = np.where(preds>0.5,'somatic','germline')
nb_germline = len(test_Y_m[test_Y_m==0])
nb_somatic = len(test_Y_m[test_Y_m==1])
nb_germline_pred = len(test_pred_m[test_pred_m['tabnet_preds']=='germline'])
nb_somatic_pred = len(test_pred_m[test_pred_m['tabnet_preds']=='somatic'])
print('[test melanoma set] real nb of germline: ' 
      + str(nb_germline) + ', somatic: ' 
      + str(nb_somatic))
print('[test melanoma set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test mixtcga
preds = classifier.predict_proba(X_test_mix)[:, 1]
test_pred_mix['tabnet_preds'] = np.where(preds>0.5,'somatic','germline')
nb_germline = len(test_Y_mix[test_Y_mix==0])
nb_somatic = len(test_Y_mix[test_Y_mix==1])
nb_germline_pred = len(test_pred_mix[test_pred_mix['tabnet_preds']=='germline'])
nb_somatic_pred = len(test_pred_mix[test_pred_mix['tabnet_preds']=='somatic'])
print('[test mixtcga set] real nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))
print('[test mixtcga set] predicted nb of germline: ' 
      + str(nb_germline_pred) 
      + ', somatic: ' 
      + str(nb_somatic_pred))

### test tnbc data
preds =  classifier.predict_proba(X_test)[:, 1]
test_pred['tabnet_preds'] = np.where(preds>0.5,'somatic','germline')
nb_germline = len(test_pred[test_pred['tabnet_preds']=='germline'])
nb_somatic = len(test_pred[test_pred['tabnet_preds']=='somatic'])
print('[test TNBC set] predicted nb of germline: ' 
      + str(nb_germline) 
      + ', somatic: ' 
      + str(nb_somatic))

### save predicts output
train_pred.to_csv(output_dir + '/preds-train.csv', index = False)
validation_pred.to_csv(output_dir + '/preds-validation.csv', index = False)
test_pred_m.to_csv(output_dir + '/preds-test-melanoma.csv', index = False)
test_pred_mix.to_csv(output_dir + '/preds-test-mixtcga.csv', index = False)
test_pred.to_csv(output_dir + '/preds-test-tnbc.csv', index = False)

