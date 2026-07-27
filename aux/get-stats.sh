#!/bin/bash

all_dirs="run-70 run-71 run-72 run-73 run-74 run-75 run-76 run-77 run-78 \
run-79 run-80 run-81 run-82 run-83 run-84 run-85 run-86 run-87 run-88 \
run-89 run-90 run-91 run-59"

cd ..

printf "N_run\tN_features\
\tauc_validation_xgboost\tauc_test_mm_xgboost\tauc_test_tcga_xgboost\
\taccuracy_validation_xgboost\tprecision_validation_xgboost\
\trecall_validation_xgboost\ttnr_validation_xgboost\
\tnpv_validation_xgboost\taccuracy_test_melanoma_xgboost\
\tprecision_test_melanoma_xgboost\trecall_test_melanoma_xgboost\
\ttnr_test_melanoma_xgboost\tnpv_test_melanoma_xgboost\
\taccuracy_test_mixtcga_xgboost\tprecision_test_mixtcga_xgboost\
\trecall_test_mixtcga_xgboost\ttnr_test_mixtcga_xgboost\
\tnpv_test_mixtcga_xgboost\tauc_validation_lightgbm\
\tauc_test_mm_lightgbm\tauc_test_tcga_lightgbm\
\taccuracy_validation_lightgbm\tprecision_validation_lightgbm\
\trecall_validation_lightgbm\ttnr_validation_lightgbm\
\tnpv_validation_lightgbm\taccuracy_test_melanoma_lightgbm\
\tprecision_test_melanoma_lightgbm\trecall_test_melanoma_lightgbm\
\ttnr_test_melanoma_lightgbm\tnpv_test_melanoma_lightgbm\
\taccuracy_test_mixtcga_lightgbm\tprecision_test_mixtcga_lightgbm\
\trecall_test_mixtcga_lightgbm\ttnr_test_mixtcga_lightgbm\
\tnpv_test_mixtcga_lightgbm\tauc_validation_tabnet\
\tauc_test_mm_tabnet\tauc_test_tcga_tabnet\taccuracy_validation_tabnet\
\tprecision_validation_tabnet\trecall_validation_tabnet\ttnr_validation_tabnet\
\tnpv_validation_tabnet\taccuracy_test_melanoma_tabnet\
\tprecision_test_melanoma_tabnet\trecall_test_melanoma_tabnet\
\ttnr_test_melanoma_tabnet\tnpv_test_melanoma_tabnet\
\taccuracy_test_mixtcga_tabnet\tprecision_test_mixtcga_tabnet\
\trecall_test_mixtcga_tabnet\ttnr_test_mixtcga_tabnet\
\tnpv_test_mixtcga_tabnet\n"

check_cols() {
	local file=${1} start=${2} end=${3} model=${4} run=${5}
	local bad
	bad=$(head -1 "${file}" | cut -f ${start}-${end} | awk -F'\t' -v m="${model}" '{for(i=1;i<=NF;i++) if ($i !~ m) print i}')
	if [ -n "${bad}" ]; then
		echo "ERROR: ${run}: columns ${start}-${end} in ${file} are not all '${model}' (mismatched field(s): ${bad})" >&2
		exit 1
	fi
}

for i in ${all_dirs}; do
	n_run=${i#run-}
	res=${i}/res
	scores=${res}/scores.txt

	check_cols "${scores}" 10 14 "XGBoost" "${i}"
	check_cols "${scores}" 21 25 "LightGBM" "${i}"
	check_cols "${scores}" 32 36 "TabNet" "${i}"

	max_feature_idx=$(grep -m1 "^max_feature_idx=" ${res}/model-lgbm.txt | cut -d= -f2)
	n_features=$((max_feature_idx + 1))

	auc_val_xgb=$(cat ${res}/auc-xgboost-validation.txt)
	auc_mm_xgb=$(cat ${res}/auc-xgboost-test-melanoma.txt)
	auc_tcga_xgb=$(cat ${res}/auc-xgboost-test-mixtcga.txt)

	auc_val_lgb=$(cat ${res}/auc-lgbm-validation.txt)
	auc_mm_lgb=$(cat ${res}/auc-lgbm-test-melanoma.txt)
	auc_tcga_lgb=$(cat ${res}/auc-lgbm-test-mixtcga.txt)

	auc_val_tab=$(cat ${res}/auc-tabnet-validation.txt)
	auc_mm_tab=$(cat ${res}/auc-tabnet-test-melanoma.txt)
	auc_tcga_tab=$(cat ${res}/auc-tabnet-test-mixtcga.txt)

	# columns 10-14: precision, recall, accuracy, tnr, npv (xgboost)
	read -r prec_val_xgb rec_val_xgb acc_val_xgb tnr_val_xgb npv_val_xgb <<< "$(cut -f 10-14 ${scores} | sed -n 3p)"
	read -r prec_mel_xgb rec_mel_xgb acc_mel_xgb tnr_mel_xgb npv_mel_xgb <<< "$(cut -f 10-14 ${scores} | sed -n 4p)"
	read -r prec_tcga_xgb rec_tcga_xgb acc_tcga_xgb tnr_tcga_xgb npv_tcga_xgb <<< "$(cut -f 10-14 ${scores} | sed -n 5p)"

	# columns 21-25: precision, recall, accuracy, tnr, npv (lightgbm)
	read -r prec_val_lgb rec_val_lgb acc_val_lgb tnr_val_lgb npv_val_lgb <<< "$(cut -f 21-25 ${scores} | sed -n 3p)"
	read -r prec_mel_lgb rec_mel_lgb acc_mel_lgb tnr_mel_lgb npv_mel_lgb <<< "$(cut -f 21-25 ${scores} | sed -n 4p)"
	read -r prec_tcga_lgb rec_tcga_lgb acc_tcga_lgb tnr_tcga_lgb npv_tcga_lgb <<< "$(cut -f 21-25 ${scores} | sed -n 5p)"

	# columns 32-36: precision, recall, accuracy, tnr, npv (tabnet)
	read -r prec_val_tab rec_val_tab acc_val_tab tnr_val_tab npv_val_tab <<< "$(cut -f 32-36 ${scores} | sed -n 3p)"
	read -r prec_mel_tab rec_mel_tab acc_mel_tab tnr_mel_tab npv_mel_tab <<< "$(cut -f 32-36 ${scores} | sed -n 4p)"
	read -r prec_tcga_tab rec_tcga_tab acc_tcga_tab tnr_tcga_tab npv_tcga_tab <<< "$(cut -f 32-36 ${scores} | sed -n 5p)"

	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t'\
'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t'\
'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t'\
'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t'\
'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t'\
'%s\t%s\t%s\t%s\t%s\t%s\n' \
		"${n_run}" "${n_features}" \
		"${auc_val_xgb}" "${auc_mm_xgb}" "${auc_tcga_xgb}" \
		"${acc_val_xgb}" "${prec_val_xgb}" "${rec_val_xgb}" \
		"${tnr_val_xgb}" "${npv_val_xgb}" \
		"${acc_mel_xgb}" "${prec_mel_xgb}" "${rec_mel_xgb}" \
		"${tnr_mel_xgb}" "${npv_mel_xgb}" \
		"${acc_tcga_xgb}" "${prec_tcga_xgb}" "${rec_tcga_xgb}" \
		"${tnr_tcga_xgb}" "${npv_tcga_xgb}" \
		"${auc_val_lgb}" "${auc_mm_lgb}" "${auc_tcga_lgb}" \
		"${acc_val_lgb}" "${prec_val_lgb}" "${rec_val_lgb}" \
		"${tnr_val_lgb}" "${npv_val_lgb}" \
		"${acc_mel_lgb}" "${prec_mel_lgb}" "${rec_mel_lgb}" \
		"${tnr_mel_lgb}" "${npv_mel_lgb}" \
		"${acc_tcga_lgb}" "${prec_tcga_lgb}" "${rec_tcga_lgb}" \
		"${tnr_tcga_lgb}" "${npv_tcga_lgb}" \
		"${auc_val_tab}" "${auc_mm_tab}" "${auc_tcga_tab}" \
		"${acc_val_tab}" "${prec_val_tab}" "${rec_val_tab}" \
		"${tnr_val_tab}" "${npv_val_tab}" \
		"${acc_mel_tab}" "${prec_mel_tab}" "${rec_mel_tab}" \
		"${tnr_mel_tab}" "${npv_mel_tab}" \
		"${acc_tcga_tab}" "${prec_tcga_tab}" "${rec_tcga_tab}" \
		"${tnr_tcga_tab}" "${npv_tcga_tab}"
done
