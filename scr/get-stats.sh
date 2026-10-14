#!/bin/bash

### A script that generates the input table for "plot-featprog-reps.R",
### with metrics: auc, accuracy, precision, recall, tnr, npv.
### It needs run-<N><replicate> folders, e.g. ../run-1a ... ../run-23e,
### next to the folder holding this script, each with the results of its run:
### ../run-1a/res/auc-*.txt
### ../run-1a/res/scores.txt
### ../run-1a/res/model-lgbm.txt
### Usage: get-stats.sh <split> [<split> ...] > stats-feat-by-feat-reps-tests.tsv
### e.g.:  get-stats.sh test-melanoma test-mixtcga test-mixanc test-ffpe-dlbcl

## settings -------------------------------------------------------------------

### every run-<N><replicate> folder (e.g. run-12c) next to the folder containing
### this script, in run order; any other run-* folder is left out
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
all_dirs=$(cd "${script_dir}/.." && for d in run-*/; do
	d=${d%/}
	[[ -d "${d}" && "${d}" =~ ^run-[0-9]+[a-e]$ ]] && printf '%s\n' "${d}"
done | sort -V | tr '\n' ' ')

# set names as in scores.txt and in the auc file names, from the command line
all_splits="${*}"

### <column label>:<auc file tag>:<first scores column>: ...
### <last scores column>:<model name in scores.txt>
### scores columns hold precision, recall, accuracy, tnr, npv
all_models="xgboost:xgboost:10:14:XGBoost \
lightgbm:lgbm:21:25:LightGBM \
tabnet:tabnet:32:36:TabNet"

## clmnt ----------------------------------------------------------------------

cd "${script_dir}/.."

tab=$'\t'

die() {
	echo "ERROR: ${1}" >&2
	exit 1
}

require_file() {
	local file=${1} run=${2}
	[ -f "${file}" ] || die "${run}: missing file ${file}"
}

check_cols() {
	local file=${1} start=${2} end=${3} model=${4} run=${5}
	local bad
	bad=$(head -1 "${file}" | cut -f ${start}-${end} | \
	  awk -F'\t' -v m="${model}" '{for(i=1;i<=NF;i++) if ($i !~ m) print i}')
	if [ -n "${bad}" ]; then
		echo "ERROR: ${run}: columns ${start}-${end} in ${file} are" \
		  "not all '${model}' (mismatched field(s): ${bad})" >&2
		exit 1
	fi
}

### a split missing from a run (e.g. runs older than test-mixanc or
### test-ffpe-dlbcl) has neither its auc files nor its row in scores.txt, so
### both readers fall back to NA
read_auc() {
	if [ -f "${1}" ]; then cat "${1}"; else echo NA; fi
}

### scores row picked by set name rather than line number
read_set() {
	local file=${1} start=${2} end=${3} set=${4}
	awk -F'\t' -v s="${set}" -v a="${start}" -v b="${end}" '
		$1 == s { for (i = a; i <= b; i++) printf "%s%s", $i, (i < b ? " " : "\n"); f = 1; exit }
		END { if (!f) print "NA NA NA NA NA" }' "${file}"
}

[ -n "${all_splits}" ] || die "no split given; usage: $(basename "${0}") <split> [<split> ...]"
[ -n "${all_dirs}" ] || die "no run-<N><a-e> folders in $(pwd)"

### a split found in no run at all is most likely a typo
for split in ${all_splits}; do
	found=
	for i in ${all_dirs}; do
		awk -F'\t' -v s="${split}" '$1 == s { f = 1; exit } END { exit !f }' \
		  "${i}/res/scores.txt" 2>/dev/null && { found=1; break; }
	done
	[ -n "${found}" ] || die "split '${split}' not found in the scores.txt of any run"
done

### header: <metric>_<split>_<model>, e.g. auc_test_melanoma_xgboost
header="N_run${tab}N_features"
for split in ${all_splits}; do
	for model in ${all_models}; do
		IFS=: read -r label tag start end name <<< "${model}"
		for metric in auc accuracy precision recall tnr npv; do
			header+="${tab}${metric}_${split//-/_}_${label}"
		done
	done
done
printf '%s\n' "${header}"

for i in ${all_dirs}; do
	n_run=${i#run-}
	res=${i}/res
	scores=${res}/scores.txt
	model_lgbm=${res}/model-lgbm.txt

	require_file "${scores}" "${i}"
	require_file "${model_lgbm}" "${i}"

	for model in ${all_models}; do
		IFS=: read -r label tag start end name <<< "${model}"
		check_cols "${scores}" ${start} ${end} "${name}" "${i}"
	done

	max_feature_idx=$(grep -m1 "^max_feature_idx=" ${model_lgbm} | cut -d= -f2)
	[[ "${max_feature_idx}" =~ ^[0-9]+$ ]] || \
	  die "${i}: missing or non-numeric 'max_feature_idx=' row in ${model_lgbm}"
	n_features=$((max_feature_idx + 1))

	row="${n_run}${tab}${n_features}"
	for split in ${all_splits}; do
		for model in ${all_models}; do
			IFS=: read -r label tag start end name <<< "${model}"
			auc=$(read_auc ${res}/auc-${tag}-${split}.txt)
			read -r prec rec acc tnr npv <<< "$(read_set ${scores} ${start} ${end} ${split})"
			row+="${tab}${auc}${tab}${acc}${tab}${prec}${tab}${rec}${tab}${tnr}${tab}${npv}"
		done
	done
	printf '%s\n' "${row}"
done
