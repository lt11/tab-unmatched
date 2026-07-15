### settings
run_tab="run-70a"
pipe_type="mulo-wesml"
dir_base="/data/ssd/lorenzot/prog"        # path to the base directory where the data is stored
cd ${dir_base}/${pipe_type}

### data

tar -czf train-vcf.tgz run-7/var-calls/gt/*.vcf.gz run-8/var-calls/gt/*.vcf.gz run-9/var-calls/gt/*.vcf.gz &
tar -czf valid-vcf.tgz run-10/var-calls/gt/*.vcf.gz run-12/var-calls/gt/*.vcf.gz &
tar -czf test-metmel.tgz run-11/var-calls/gt/*.vcf.gz &
tar -czf test-mixtcga.tgz run-13/var-calls/gt/*.vcf.gz run-15/var-calls/gt/*.vcf.gz &
### different extension (mulo-wesunc)
tar -czf test-tnbc.tgz ../tnbc/run-1/var-calls/anno-snpeff-unc/*dbs.vcf.gz &

wait
mv *tgz ${dir_base}/tab-ml

### tab-ml

cd ${dir_base}/tab-ml
mkdir -p ${dir_base}/tab-ml/${run_tab}/train
mkdir -p ${dir_base}/tab-ml/${run_tab}/validation
mkdir -p ${dir_base}/tab-ml/${run_tab}/test/melanoma
mkdir -p ${dir_base}/tab-ml/${run_tab}/test/mixtcga
mkdir -p ${dir_base}/tab-ml/${run_tab}/test/tnbc

cd ${dir_base}/tab-ml/
mv train-vcf.tgz ${dir_base}/tab-ml/${run_tab}/train
cd ${dir_base}/tab-ml/${run_tab}/train
tar -zvxf train-vcf.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf train-vcf.tgz 
gunzip *vcf.gz

cd ${dir_base}/tab-ml/
mv valid-vcf.tgz ${dir_base}/tab-ml/${run_tab}/validation
cd ${dir_base}/tab-ml/${run_tab}/validation
tar -zvxf valid-vcf.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf valid-vcf.tgz
gunzip *vcf.gz

cd ${dir_base}/tab-ml/
mv test-metmel.tgz ${dir_base}/tab-ml/${run_tab}/test/melanoma
cd ${dir_base}/tab-ml/${run_tab}/test/melanoma
tar -zvxf test-metmel.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf test-metmel.tgz
gunzip *vcf.gz

cd ${dir_base}/tab-ml/
mv test-mixtcga.tgz ${dir_base}/tab-ml/${run_tab}/test/mixtcga
cd ${dir_base}/tab-ml/${run_tab}/test/mixtcga
tar -zvxf test-mixtcga.tgz  
mv run-*/var-calls/gt/*vcf.gz .
rm -rf test-mixtcga.tgz
gunzip *vcf.gz

cd ${dir_base}/tab-ml
mv test-tnbc.tgz ${dir_base}/tab-ml/${run_tab}/test/tnbc
cd ${dir_base}/tab-ml/${run_tab}/test/tnbc
tar -zvxf test-tnbc.tgz
mv tnbc/run-*/var-calls/anno-snpeff-unc/*dbs.vcf.gz .
rm -rf test-tnbc.tgz
gunzip *vcf.gz

cd ${dir_base}  
dir_tab=$(find ./tab-ml -name tab | tail -1)
cp -r ${dir_tab} ./tab-ml/${run_tab}
cd ./tab-ml/${run_tab}

### set seed
all_seeds=(21709 11796 17856 24207 12554)

### run a
ind_run=0
cd scr
python unc-normal-drop-5.py ${all_seeds[ind_run]} > log-unc-normal-drop-5.txt
cd ../aux
Rscript make-score.R > log-make-score.txt
Rscript calculate-tmb.R 
Rscript box-tmb.R 
Rscript lin-regression-tmb.R
Rscript make-score-snv.R
Rscript make-score-sid.R
Rscript make-score-sbs.R
Rscript make-score-tbt.R

### replicate runs
run_pref="${run_tab%?}"

for ind_rep in b c d e; do
  dir_runs=$(cd ../.. && pwd)

  ### define the new run name, e.g. run-68b
  dir_run_new="${dir_runs}/${run_pref}${ind_rep}"

  ### find a run-Nx folder containing all required inputs
  dir_run_source=""

  for dir_candidate in "${dir_runs}/${run_pref}"?; do
    if [[ "${dir_candidate}" != "${dir_run_new}" && \
          -d "${dir_candidate}/scr" && \
          -d "${dir_candidate}/aux" && \
          -d "${dir_candidate}/train" && \
          -d "${dir_candidate}/validation" && \
          -d "${dir_candidate}/test" && \
          -d "${dir_candidate}/tab" ]]; then
      dir_run_source="${dir_candidate}"
      break
    fi
  done

  if [[ -z "${dir_run_source}" ]]; then
    echo "Error: no ${run_pref}x folder contains input data." >&2
    exit 1
  fi

  printf 'Preparing %s from %s\n' \
    "$(basename "${dir_run_new}")" \
    "$(basename "${dir_run_source}")"

  ### create the new run folder
  mkdir -p "${dir_run_new}/scr" "${dir_run_new}/aux"

  ### copy code into the new run folder
  cp -a "${dir_run_source}/scr/." "${dir_run_new}/scr/"
  cp -a "${dir_run_source}/aux/." "${dir_run_new}/aux/"

  ### move variants and table into the new run folder
  mv "${dir_run_source}/train" "${dir_run_new}/"
  mv "${dir_run_source}/validation" "${dir_run_new}/"
  mv "${dir_run_source}/test" "${dir_run_new}/"
  mv "${dir_run_source}/tab" "${dir_run_new}/"

  ### set seed and run
  ind_run=$((ind_run + 1))
  echo "Running $(basename "${dir_run_new}") with seed ${all_seeds[ind_run]}"
  cd "${dir_run_new}/scr"
  python unc-normal-drop-5.py ${all_seeds[ind_run]} > log-unc-normal-drop-5.txt

  cd "${dir_run_new}/aux"
  Rscript make-score.R > log-make-score.txt
  Rscript calculate-tmb.R
  Rscript box-tmb.R
  Rscript lin-regression-tmb.R
  Rscript make-score-snv.R
  Rscript make-score-sid.R
  Rscript make-score-sbs.R
  Rscript make-score-tbt.R
done
