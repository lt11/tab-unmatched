#!/bin/bash

### Runs a whole feature progression: run-1<r> keeps every feature, run-2<r>
### drops the first feature of drop_order, run-3<r> the first two and so on,
### with one seeded replicate <r> per letter in reps; then it collects the
### stats with get-stats.sh and plots them with plot-featprog-reps.R.
### It sits in ${dir_base}/scr next to the rest of the code; metadata.csv and
### index-paths.txt, both kept by hand, go in ${dir_base}/dat/tab.
### Usage: nohup bash scr/superrunner.sh > superrunner.log 2>&1 &
### It stops at the first failure; a relaunch skips what is done already
### (the data, the finished runs) and redoes an unfinished run from scratch.

## settings -------------------------------------------------------------------

### where the mulo-wesml VCFs are, the first of these folders that exists;
### the paths in index-paths.txt are relative to it
vcf_src_candidates=(/data/ssd/lorenzot/prog/mulo-wesml
                    /Users/Lorenzo/prog/indi-cange/mulo-wesml)
### logistic regression and hard filters
run_logit="true"
run_hf="true"
### the replicates and their seeds, in the same order
reps=(a b c d e)
all_seeds=(21709 11796 17856 24207 12554)
### the features dropped cumulatively, run-N drops the first N-1: any of
### valid_drop_features in the python script, quoted if holding '>'
drop_order=(CNT
            ontology_NR
            DP
            ontology_inframe_indel
            AF
            ontology_frame_shift_indel
            AD
            pop_max
            subs_type_NA
            ontology_nonsense
            'subs_type_G>A'
            trinucleotide_TCA
            trinucleotide_GGG
            'subs_type_G>T'
            'subs_type_C>T'
            ontology_missense
            'subs_type_G>C'
            trinucleotide_TCC
            trinucleotide_GGA
            trinucleotide_TGA
            trinucleotide_GTG
            'subs_type_C>A')
### the interpreters
py_bin="python"
r_bin="Rscript"

### fixed settings: the code, all of it next to this script
py_main="unc-normal-drop-9.py"
r_per_run=(make-score calculate-tmb box-tmb
           make-score-snv make-score-sid make-score-sbs make-score-tbt)
### the sets of the split column, each with its folder under dat/vcf, in the
### order they are reported
sets_all="train:train validation:validation test-mixtcga:test/mixtcga \
test-mixanc:test/mixanc test-metmel:test/melanoma \
test-ffpe-dlbcl:test/ffpe-dlbcl"
### the same sets as the python script names them in its results; it writes
### an auc file per model for every set but train
sets_py="train validation test-melanoma test-mixtcga test-ffpe-dlbcl \
test-mixanc"
auc_models="xgboost lgbm tabnet"
### the splits handed to get-stats.sh
splits_test="test-melanoma test-mixtcga test-mixanc test-ffpe-dlbcl"
splits_trval="train validation"

### folders: the experiment is the folder above the one holding this script
dir_scr=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
dir_base=$(dirname "${dir_scr}")
dir_dat="${dir_base}/dat"
file_meta="${dir_dat}/tab/metadata.csv"
file_index="${dir_dat}/tab/index-paths.txt"
dir_lst="${dir_dat}/splits"
dir_vcf="${dir_dat}/vcf"
dir_feat="${dir_dat}/feat"
dir_res="${dir_base}/res"
dir_lock="${dir_base}/.superrunner.lock"
### the run folders, e.g. run-12c
re_run='^run-([0-9]+)([a-e])$'

## functions ------------------------------------------------------------------

say () {
  printf '[superrunner] %s %s\n' "$(date '+%F %T')" "$*"
}

die () {
  printf '[superrunner] %s Error: %s\n' "$(date '+%F %T')" "$*" >&2
  exit 1
}

### checksum and size of a file, as one word
sum_of () {
  cksum < "$1" | awk '{ print $1 "-" $2 }'
}

### rm -rf fenced to what this script rebuilds: a run folder of this
### experiment, or one of the temporary data folders
clear_path () {
  if [ "${1%/*}" = "${dir_base}" ] && [[ "${1##*/}" =~ ${re_run} ]]; then
    rm -rf "$1"
  elif [ "$1" = "${dir_dat}/vcf.tmp" ] || [ "$1" = "${dir_dat}/feat.tmp" ]; then
    rm -rf "$1"
  else
    echo "Error: refusing to clear $1" >&2
    return 1
  fi
}

### copy one set straight from its path list into its folder and unpack it
copy_set () {
  local label="$1" dest="$2" lst="${dir_lst}/$1-paths.txt"
  local expect n_vcf path
  expect=$(grep -c . "${lst}" || true)
  mkdir -p "${dest}" || return 1
  while IFS= read -r path; do
    cp "${dir_vcf_src}/${path}" "${dest}/" || return 1
  done < "${lst}"
  if ls "${dest}"/*.vcf.gz > /dev/null 2>&1; then
    gunzip -f "${dest}"/*.vcf.gz || return 1
  fi
  n_vcf=$(find "${dest}" -maxdepth 1 -name '*.vcf' | wc -l | tr -d ' ')
  printf '%-16s %3d vcf in %s\n' "${label}" "${n_vcf}" "${dest#${dir_dat}/}"
  ### the folder is new, so a mismatch means the copy itself is short of
  ### what the split asked for
  if [ "${n_vcf}" -ne "${expect}" ]; then
    echo "Error: expected ${expect} vcf in ${dest}, found ${n_vcf}" >&2
    return 1
  fi
}

### the split column of metadata.csv assigns every tumour sample to one set;
### VCF file names are the tumour fastq id and normals have no VCF of their
### own; the sets are copied into dat/vcf.tmp, renamed dat/vcf once complete
prep_vcf () {
  local d pair label path bad_index n_want n_got n_gone n_dup n_abs
  local n_rep n_unknown tot_missing=0 tot_dup=0 tot_absent=0 tot_empty=0
  local dir_tmp="${dir_dat}/vcf.tmp"

  ### the source of the VCFs is needed only here
  dir_vcf_src=""
  for d in "${vcf_src_candidates[@]}"; do
    if [ -d "${d}" ]; then
      dir_vcf_src="${d}"
      break
    fi
  done
  [ -n "${dir_vcf_src}" ] \
    || die "none of the VCF sources exists:" "${vcf_src_candidates[@]}"
  say "building ${dir_vcf} from ${dir_vcf_src}"
  mkdir -p "${dir_lst}" || die "could not create ${dir_lst}"

  ### every non-blank line of the index must be a path to a .vcf.gz relative
  ### to the source, so that the same index serves every machine
  bad_index=$(awk '{ sub(/[ \t\r]+$/, "") } NF && (/^\// || !/\.vcf\.gz$/)' \
    "${file_index}")
  if [ -n "${bad_index}" ]; then
    echo "Error: not a relative path to a .vcf.gz in ${file_index}:" >&2
    printf '%s\n' "${bad_index}" >&2
    exit 1
  fi

  ### one fastq id list per set, straight from the split column; the columns
  ### are located by name and both sit ahead of the only quoted fields in the
  ### table, so splitting on commas is safe here
  for label in ${set_labels}; do
    : > "${dir_lst}/${label}-ids.txt"
  done

  awk -F',' -v dir="${dir_lst}" -v sets="${set_labels}" '
    BEGIN { n_set = split(sets, lab, " ")
            for (i = 1; i <= n_set; i++) want[lab[i]] = 1 }
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        gsub(/^[ \t]+|[ \t\r]+$/, "", $i)
        if ($i == "fastq id") c_id = i
        if ($i == "split") c_set = i
      }
      if (!c_id || !c_set) {
        print "Error: no fastq id or split column in the metadata" > "/dev/stderr"
        exit 1
      }
      next
    }
    {
      id = $c_id; set = $c_set
      gsub(/^[ \t]+|[ \t\r]+$/, "", id)
      gsub(/^[ \t]+|[ \t\r]+$/, "", set)
      if (set in want) { print id > (dir "/" set "-ids.txt"); n[set]++ }
    }
    END { for (i = 1; i <= n_set; i++)
            printf "%-16s %3d samples\n", lab[i], n[lab[i]] + 0 }
  ' "${file_meta}" || exit 1

  ### resolve the ids to VCF paths, then check that nothing is missing, double
  ### or listed in the index without being on disk; the first file of each
  ### awk is told apart by name, as an empty one would fool NR == FNR
  for label in ${set_labels}; do
    awk -F'/' 'FILENAME == ARGV[1] { want[$0] = 1; next }
               { sub(/[ \t\r]+$/, "") }
               NF { id = $NF; sub(/\.vcf\.gz$/, "", id); if (id in want) print }' \
      "${dir_lst}/${label}-ids.txt" "${file_index}" \
      > "${dir_lst}/${label}-paths.txt"
    awk -F'/' 'FILENAME == ARGV[1] { id = $NF; sub(/\.vcf\.gz$/, "", id)
                                     have[id] = 1; next }
               !($0 in have)' \
      "${dir_lst}/${label}-paths.txt" "${dir_lst}/${label}-ids.txt" \
      > "${dir_lst}/${label}-missing.txt"
    awk -F'/' '{ id = $NF; sub(/\.vcf\.gz$/, "", id); print id }' \
      "${dir_lst}/${label}-paths.txt" | sort | uniq -d \
      > "${dir_lst}/${label}-duplicated.txt"
    while IFS= read -r path; do
      [ -f "${dir_vcf_src}/${path}" ] || printf '%s\n' "${path}"
    done < "${dir_lst}/${label}-paths.txt" > "${dir_lst}/${label}-absent.txt"
    n_want=$(grep -c . "${dir_lst}/${label}-ids.txt" || true)
    n_got=$(grep -c . "${dir_lst}/${label}-paths.txt" || true)
    n_gone=$(grep -c . "${dir_lst}/${label}-missing.txt" || true)
    n_dup=$(grep -c . "${dir_lst}/${label}-duplicated.txt" || true)
    n_abs=$(grep -c . "${dir_lst}/${label}-absent.txt" || true)
    printf '%-16s %3d requested, %3d found, %3d missing, %3d duplicated, %3d absent\n' \
      "${label}" "${n_want}" "${n_got}" "${n_gone}" "${n_dup}" "${n_abs}"
    tot_missing=$((tot_missing + n_gone))
    tot_dup=$((tot_dup + n_dup))
    tot_absent=$((tot_absent + n_abs))
    ### an empty set means the split column does not use this label at all
    if [ "${n_want}" -eq 0 ]; then
      echo "Error: no sample carries split ${label} in ${file_meta}" >&2
      tot_empty=$((tot_empty + 1))
    fi
  done

  ### a fastq id on two rows of the split column, in one set or in two,
  ### would be counted twice, or trained and tested on at once
  for label in ${set_labels}; do
    cat "${dir_lst}/${label}-ids.txt"
  done | sort | uniq -d > "${dir_lst}/ids-repeated.txt"
  n_rep=$(grep -c . "${dir_lst}/ids-repeated.txt" || true)

  ### an index entry that no row of the metadata knows about ends up in no set;
  ### harmless, but most likely a typo in one of the two files
  awk 'FILENAME == ARGV[1] { sub(/[ \t\r]+$/, "")
                             if (NF) { n = split($0, p, "/"); id = p[n]
                                       sub(/\.vcf\.gz$/, "", id); idx[id] = 1 }
                             next }
       { n = split($0, f, ",") }
       FNR == 1 { for (i = 1; i <= n; i++) {
                    gsub(/^[ \t]+|[ \t\r]+$/, "", f[i])
                    if (f[i] == "fastq id") c_id = i
                  }
                  next }
       { id = f[c_id]; gsub(/^[ \t]+|[ \t\r]+$/, "", id); known[id] = 1 }
       END { for (id in idx) if (!(id in known)) print id }' \
    "${file_index}" "${file_meta}" | sort > "${dir_lst}/index-unknown.txt"
  n_unknown=$(grep -c . "${dir_lst}/index-unknown.txt" || true)
  if [ "${n_unknown}" -gt 0 ]; then
    echo "WARNING: ${n_unknown} VCFs in ${file_index} match no fastq id in" \
         "${file_meta}; see ${dir_lst}/index-unknown.txt" >&2
  fi

  ### all five are fatal and there is no override: an incomplete set would
  ### train or test on fewer samples than the split asks for, an ambiguous
  ### one would keep whichever of its VCFs happened to be copied last, and a
  ### repeated id would count a sample twice, all of them silently; every
  ### problem is reported before exiting with error
  if [ "${tot_missing}" -gt 0 ]; then
    echo "Error: ${tot_missing} samples have no VCF in ${file_index};" \
         "see ${dir_lst}/*-missing.txt" >&2
  fi
  if [ "${tot_dup}" -gt 0 ]; then
    echo "Error: ${tot_dup} samples listed more than once in ${file_index};" \
         "see ${dir_lst}/*-duplicated.txt" >&2
  fi
  if [ "${tot_absent}" -gt 0 ]; then
    echo "Error: ${tot_absent} VCFs in ${file_index} do not exist under" \
         "${dir_vcf_src}; see ${dir_lst}/*-absent.txt" >&2
  fi
  if [ "${n_rep}" -gt 0 ]; then
    echo "Error: ${n_rep} fastq ids on more than one row of ${file_meta};" \
         "see ${dir_lst}/ids-repeated.txt" >&2
  fi
  if [ "${tot_missing}" -gt 0 ] || [ "${tot_dup}" -gt 0 ] \
     || [ "${tot_absent}" -gt 0 ] || [ "${tot_empty}" -gt 0 ] \
     || [ "${n_rep}" -gt 0 ]; then
    die "the split lists in ${dir_lst} failed the checks above"
  fi

  ### copy into a pristine temporary folder, renamed once complete, so that
  ### an interrupted copy is never taken for a finished one
  clear_path "${dir_tmp}" || exit 1
  for pair in ${sets_all}; do
    copy_set "${pair%%:*}" "${dir_tmp}/${pair#*:}" \
      || die "copying the set ${pair%%:*} failed"
  done
  printf '%s\n' "${stamp_vcf}" > "${dir_tmp}/.stamp"
  mv "${dir_tmp}" "${dir_vcf}" || die "could not rename ${dir_tmp} to ${dir_vcf}"
  say "dat/vcf built"
}

### the feature tables depend on neither the seed nor the dropped features,
### so they are built once, in dat/feat.tmp renamed dat/feat once complete,
### and every run only trains on them
prep_feat () {
  local dir_tmp="${dir_dat}/feat.tmp" s
  clear_path "${dir_tmp}" || exit 1
  mkdir -p "${dir_tmp}" || die "could not create ${dir_tmp}"
  say "building ${dir_feat}"
  ( cd "${dir_tmp}" && "${py_bin}" -u "${dir_scr}/${py_main}" \
      "${all_seeds[0]}" --stage prep --dir-vcf "${dir_vcf}" \
      --file-meta "${file_meta}" --dir-feat "${dir_tmp}" \
      > "${dir_tmp}/log-prep.txt" 2>&1 ) \
    || die "building the feature tables failed, see ${dir_tmp}/log-prep.txt"
  for s in ${sets_py}; do
    [ -s "${dir_tmp}/${s}-df-flt.csv" ] \
      || die "no ${s}-df-flt.csv in ${dir_tmp}, see ${dir_tmp}/log-prep.txt"
  done
  printf 'data=%s\ncode=%s\n' "${stamp_vcf}" "${sum_py}" > "${dir_tmp}/.stamp"
  mv "${dir_tmp}" "${dir_feat}" \
    || die "could not rename ${dir_tmp} to ${dir_feat}"
  say "dat/feat built"
}

### a checksum over the python script and the R scripts of every run
sum_of_code () {
  local r
  {
    cksum < "${dir_scr}/${py_main}"
    for r in "${r_per_run[@]}"; do
      cksum < "${dir_scr}/${r}.R"
    done
  } | cksum | awk '{ print $1 "-" $2 }'
}

### the settings a finished run must share with the current ones to be
### kept: $1 is the seed, the rest the dropped features
run_config () {
  local seed="$1"
  shift
  printf 'seed=%s\nrun_logit=%s\nrun_hf=%s\ndrop_features=%s\ndata=%s\n' \
    "${seed}" "${run_logit}" "${run_hf}" "$*" "${stamp_vcf}"
}

### the number of features a run trained on, from its LightGBM model
n_features_of () {
  awk -F'=' '$1 == "max_feature_idx" && $2 ~ /^[0-9]+$/ { print $2 + 1; f = 1
                                                          exit }
             END { exit !f }' "$1/res/model-lgbm.txt" 2> /dev/null
}

### every run trains on one feature fewer than the run before it; the first
### run met, run-1 of the first replicate, sets the count for all of them
check_n_features () {
  local name="$1" n_set="$2" dir_run="$3" n_feat
  n_feat=$(n_features_of "${dir_run}") \
    || die "${name}: no max_feature_idx in res/model-lgbm.txt"
  if [ -z "${n_all}" ]; then
    n_all=$((n_feat + n_set - 1))
  fi
  [ "${n_feat}" -eq $((n_all - n_set + 1)) ] \
    || die "${name} trained on ${n_feat} features, expected" \
           "$((n_all - n_set + 1)): ${n_all} minus $((n_set - 1)) dropped"
}

### one replicate of one feature set: $1 is the run number, $2 the replicate,
### $3 its seed, the rest the dropped features (none for run-1)
run_one () {
  local n_set="$1" rep="$2" seed="$3"
  shift 3
  local name="run-${n_set}${rep}"
  local dir_run="${dir_base}/run-${n_set}${rep}"
  local log_py="${dir_run}/log/log-${py_main%.py}.txt"
  local want have f m r s t_start r_args
  want=$(run_config "${seed}" "$@")

  ### a finished run is kept, provided it was made with the current settings
  if [ -f "${dir_run}/.done" ]; then
    have=$(grep -v '^code=' "${dir_run}/run-config.txt" 2> /dev/null)
    [ "${have}" = "${want}" ] \
      || die "${name} was made with other settings, see" \
             "${dir_run}/run-config.txt; change the settings back, or" \
             "remove ${name} to redo it"
    grep -qxF "code=${sum_code}" "${dir_run}/run-config.txt" \
      || say "WARNING: ${name} was made with other versions of the scripts"
    check_n_features "${name}" "${n_set}" "${dir_run}"
    say "${name} done already, kept"
    n_kept=$((n_kept + 1))
    return 0
  fi

  ### anything else left by an earlier launch is redone from scratch
  if [ -e "${dir_run}" ]; then
    say "${name} was left unfinished, redoing it"
  fi
  clear_path "${dir_run}" || exit 1
  mkdir -p "${dir_run}/log" || die "could not create ${dir_run}/log"
  { printf '%s\n' "${want}"; printf 'code=%s\n' "${sum_code}"; } \
    > "${dir_run}/run-config.txt"
  t_start=${SECONDS}
  say "${name}: seed ${seed}, $((n_set - 1)) features dropped"

  ### the working directory is the log folder, so that nothing the scripts
  ### might write next to it lands outside the run
  ( cd "${dir_run}/log" && "${py_bin}" -u "${dir_scr}/${py_main}" \
      "${seed}" "${run_logit}" "${run_hf}" --drop-features "$@" \
      --stage train --dir-feat "${dir_feat}" --dir-out "${dir_run}/res" \
      > "${log_py}" 2>&1 ) \
    || die "${name}: ${py_main} failed, see ${log_py}"
  for s in ${sets_py}; do
    [ -s "${dir_run}/res/preds-${s}.csv" ] \
      || die "${name}: ${py_main} wrote no res/preds-${s}.csv"
    [ "${s}" = "train" ] && continue
    for m in ${auc_models}; do
      [ -s "${dir_run}/res/auc-${m}-${s}.txt" ] \
        || die "${name}: ${py_main} wrote no res/auc-${m}-${s}.txt"
    done
  done
  check_n_features "${name}" "${n_set}" "${dir_run}"

  for r in "${r_per_run[@]}"; do
    r_args=("${dir_run}")
    if [ "${r}" = "calculate-tmb" ]; then
      r_args=("${dir_run}" "${file_meta}")
    fi
    ( cd "${dir_run}/log" && "${r_bin}" "${dir_scr}/${r}.R" "${r_args[@]}" \
        > "${dir_run}/log/log-${r}.txt" 2>&1 ) \
      || die "${name}: ${r}.R failed, see ${dir_run}/log/log-${r}.txt"
  done
  for f in scores.txt metadata-tmb.txt; do
    [ -s "${dir_run}/res/${f}" ] \
      || die "${name}: no res/${f} after the R scripts"
  done

  date '+%F %T' > "${dir_run}/.done"
  n_trained=$((n_trained + 1))
  say "${name} finished in $(((SECONDS - t_start) / 60)) min"
}

## preflight ------------------------------------------------------------------

### the experiment folder must be laid out as described at the top, which
### also stops a launch from the repository by mistake
if [ ! -f "${file_meta}" ]; then
  die "metadata.csv not found at ${file_meta}; is ${dir_base} an" \
      "experiment folder?"
fi
if ! grep -q . "${file_index}" 2> /dev/null; then
  die "index-paths.txt missing or empty at ${file_index}"
fi

### one superrunner at a time per experiment folder
if ! mkdir "${dir_lock}" 2> /dev/null; then
  die "${dir_lock} exists: another superrunner (pid" \
      "$(cat "${dir_lock}/pid" 2> /dev/null)) is running here; if none" \
      "is, remove the lock with: rm -r ${dir_lock}"
fi
trap 'rm -rf "${dir_lock}"' EXIT
trap 'say "interrupted; a relaunch redoes any run left unfinished"; exit 130' INT
trap 'say "terminated; a relaunch redoes any run left unfinished"; exit 143' TERM
echo "$$" > "${dir_lock}/pid"

### every plot goes to a file: without a display-free backend, each
### plt.show() of the python script would wait on a window wherever a
### desktop or a forwarded X display is at hand
export MPLBACKEND=Agg

for prog in "${py_bin}" "${r_bin}" gunzip cksum; do
  command -v "${prog}" > /dev/null || die "${prog} not found"
done
for f in "${py_main}" get-stats.sh plot-featprog-reps.R; do
  [ -f "${dir_scr}/${f}" ] || die "${f} not found in ${dir_scr}"
done
for r in "${r_per_run[@]}"; do
  [ -f "${dir_scr}/${r}.R" ] || die "${r}.R not found in ${dir_scr}"
done

### the replicates: letters a to e, the ones get-stats.sh reads, one seed each
[ "${#reps[@]}" -gt 0 ] || die "no replicates in reps"
[ "${#all_seeds[@]}" -eq "${#reps[@]}" ] \
  || die "${#reps[@]} replicates but ${#all_seeds[@]} seeds"
for r in "${reps[@]}"; do
  case "${r}" in
    [a-e]) ;;
    *) die "replicate '${r}' is not one of a b c d e" ;;
  esac
done
[ -z "$(printf '%s\n' "${reps[@]}" | sort | uniq -d)" ] \
  || die "a replicate is listed twice in reps"
for s in "${all_seeds[@]}"; do
  [[ "${s}" =~ ^[0-9]+$ ]] || die "seed '${s}' is not a number"
done

### the features to drop: no repeats, no 'none', and every name known to the
### python script, which refuses the unknown ones
for f in "${drop_order[@]}"; do
  [ "${f}" != "none" ] \
    || die "drop_order holds 'none', but run-1 keeps every feature already"
done
dup_drops=$(printf '%s\n' "${drop_order[@]}" | sort | uniq -d)
[ -z "${dup_drops}" ] \
  || die "listed more than once in drop_order:" ${dup_drops}
if [ "${#drop_order[@]}" -gt 0 ]; then
  "${py_bin}" "${dir_scr}/${py_main}" "${all_seeds[0]}" --dry-run \
    --drop-features "${drop_order[@]}" > /dev/null \
    || die "${py_main} refused drop_order, see above"
fi

### a run folder outside the current plan would end up in the stats
n_sets=$((${#drop_order[@]} + 1))
for d in "${dir_base}"/run-*; do
  name="${d##*/}"
  [ -d "${d}" ] && [[ "${name}" =~ ${re_run} ]] || continue
  n="${BASH_REMATCH[1]}"
  r="${BASH_REMATCH[2]}"
  known=""
  for x in "${reps[@]}"; do
    if [ "${x}" = "${r}" ]; then
      known="yes"
    fi
  done
  if [ "${n}" -lt 1 ] || [ "${n}" -gt "${n_sets}" ] || [ -z "${known}" ]; then
    die "${name} is not part of the current plan (run-1 to run-${n_sets}," \
        "replicates ${reps[*]}); remove it, or change the settings back"
  fi
done

say "${dir_base}: ${n_sets} feature sets x ${#reps[@]} replicates" \
    "= $((n_sets * ${#reps[@]})) runs"

## variants -------------------------------------------------------------------

set_labels=""
for pair in ${sets_all}; do
  set_labels="${set_labels:+${set_labels} }${pair%%:*}"
done
stamp_vcf="meta=$(sum_of "${file_meta}") index=$(sum_of "${file_index}")"

if [ -d "${dir_vcf}" ]; then
  [ -f "${dir_vcf}/.stamp" ] \
    || die "${dir_vcf} has no .stamp, so this script did not complete it;" \
           "remove it to rebuild it"
  [ "$(cat "${dir_vcf}/.stamp")" = "${stamp_vcf}" ] \
    || die "${dir_vcf} was built from another metadata.csv or" \
           "index-paths.txt; to rebuild it, remove ${dir_vcf}, ${dir_feat}" \
           "and every run-* folder"
  say "dat/vcf is up to date"
else
  prep_vcf
fi

## features -------------------------------------------------------------------

sum_py=$(sum_of "${dir_scr}/${py_main}")
if [ -d "${dir_feat}" ]; then
  [ -f "${dir_feat}/.stamp" ] \
    || die "${dir_feat} has no .stamp, so this script did not complete it;" \
           "remove it to rebuild it"
  grep -qxF "data=${stamp_vcf}" "${dir_feat}/.stamp" \
    || die "${dir_feat} was built from other VCFs than those in ${dir_vcf};" \
           "to rebuild it, remove ${dir_feat} and every run-* folder"
  grep -qxF "code=${sum_py}" "${dir_feat}/.stamp" \
    || say "WARNING: ${dir_feat} was built by another version of ${py_main}"
  say "dat/feat is up to date"
else
  prep_feat
fi

### every feature to drop must be a column of the training table, or its
### run would train on one feature more than planned
header=$(head -1 "${dir_feat}/train-df-flt.csv" | tr -d '"\r')
for f in "${drop_order[@]}"; do
  case ",${header}," in
    *",${f},"*) ;;
    *) die "${f} in drop_order is not a column of" \
           "${dir_feat}/train-df-flt.csv" ;;
  esac
done

## clmnt ----------------------------------------------------------------------

sum_code=$(sum_of_code)
n_all=""
n_trained=0
n_kept=0

n_set=1
while [ "${n_set}" -le "${n_sets}" ]; do
  ### run-N drops the first N-1 features of drop_order
  if [ "${n_set}" -eq 1 ]; then
    drops=(none)
  else
    drops=("${drop_order[@]:0:$((n_set - 1))}")
  fi
  i_rep=0
  while [ "${i_rep}" -lt "${#reps[@]}" ]; do
    run_one "${n_set}" "${reps[i_rep]}" "${all_seeds[i_rep]}" "${drops[@]}"
    i_rep=$((i_rep + 1))
  done
  n_set=$((n_set + 1))
done

## stats ----------------------------------------------------------------------

mkdir -p "${dir_res}/accu-prog" "${dir_res}/accu-prog-plots" \
  || die "could not create the result folders in ${dir_res}"
file_test="${dir_res}/accu-prog/stats-feat-by-feat-reps-tests.tsv"
file_trval="${dir_res}/accu-prog/stats-feat-by-feat-reps-train-validation.tsv"
file_log_plot="${dir_res}/log-plot-featprog-reps.txt"

### each table is written aside and renamed, so that a failure never leaves
### a partial one behind
say "stats on ${splits_test}"
bash "${dir_scr}/get-stats.sh" ${splits_test} > "${file_test}.tmp" \
  || die "get-stats.sh failed on ${splits_test}"
mv "${file_test}.tmp" "${file_test}" || die "could not write ${file_test}"

say "plots of the test stats"
"${r_bin}" "${dir_scr}/plot-featprog-reps.R" "${file_test}" \
  "${dir_res}/accu-prog-plots" > "${file_log_plot}" 2>&1 \
  || die "plot-featprog-reps.R failed, see ${file_log_plot}"

say "stats on ${splits_trval}"
bash "${dir_scr}/get-stats.sh" ${splits_trval} > "${file_trval}.tmp" \
  || die "get-stats.sh failed on ${splits_trval}"
mv "${file_trval}.tmp" "${file_trval}" || die "could not write ${file_trval}"

say "all done: ${n_trained} runs trained, ${n_kept} kept from before, in" \
    "$((SECONDS / 3600)) h $((SECONDS % 3600 / 60)) min"
say "stats in ${dir_res}/accu-prog, plots in ${dir_res}/accu-prog-plots"
