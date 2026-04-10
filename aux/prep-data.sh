### settings
run_tab="run-68"
pipe_type="mulo-wesone"
${dir_base}=EDIT-ME        # path to the base directory where the data is stored
cd ${dir_base}/$pipe_type

### data
tar -czf train-vcf.tgz run-7/var-calls/gt/*.vcf.gz run-8/var-calls/gt/*.vcf.gz run-9/var-calls/gt/*.vcf.gz &
tar -czf valid-vcf.tgz run-10/var-calls/gt/*.vcf.gz run-12/var-calls/gt/*.vcf.gz &
tar -czf test-metmel.tgz run-11/var-calls/gt/*.vcf.gz &
tar -czf test-mixtcga.tgz run-13/var-calls/gt/*.vcf.gz run-15/var-calls/gt/*.vcf.gz &
wait
mv *tgz ${dir_base}/tab-ml

### tab-ml

cd ${dir_base}/tab-ml
mkdir -p ${dir_base}/tab-ml/$run_tab/train
mkdir -p ${dir_base}/tab-ml/$run_tab/validation
mkdir -p ${dir_base}/tab-ml/$run_tab/test/melanoma
mkdir -p ${dir_base}/tab-ml/$run_tab/test/mixtcga
mkdir -p ${dir_base}/tab-ml/$run_tab/test/tnbc

cd ${dir_base}/tab-ml/
mv train-vcf.tgz ${dir_base}/tab-ml/$run_tab/train
cd ${dir_base}/tab-ml/$run_tab/train
tar -zvxf train-vcf.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf train-vcf.tgz 
gunzip *vcf.gz

cd ${dir_base}/tab-ml/
mv valid-vcf.tgz ${dir_base}/tab-ml/$run_tab/validation
cd ${dir_base}/tab-ml/$run_tab/validation
tar -zvxf valid-vcf.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf valid-vcf.tgz
gunzip *vcf.gz

cd ${dir_base}/tab-ml/
mv test-metmel.tgz ${dir_base}/tab-ml/$run_tab/test/melanoma
cd ${dir_base}/tab-ml/$run_tab/test/melanoma
tar -zvxf test-metmel.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf test-metmel.tgz
gunzip *vcf.gz

cd ${dir_base}/tab-ml/
mv test-mixtcga.tgz ${dir_base}/tab-ml/$run_tab/test/mixtcga
cd ${dir_base}/tab-ml/$run_tab/test/mixtcga
tar -zvxf test-mixtcga.tgz  
mv run-*/var-calls/gt/*vcf.gz .
rm -rf test-mixtcga.tgz
gunzip *vcf.gz

cd ${dir_base}/tab-ml
mv test-tnbc.tgz ${dir_base}/tab-ml/$run_tab/test/tnbc
cd ${dir_base}/tab-ml/$run_tab/test/tnbc
tar -zvxf test-tnbc.tgz
mv tnbc/run-*/var-calls/anno-snpeff-unc/*dbs.vcf.gz .
rm -rf test-tnbc.tgz
gunzip *vcf.gz

cd ${dir_base}  
dir_tab=$(find ./tab-ml -name tab | tail -1)
cp -r ${dir_tab} ./tab-ml/$run_tab
cd ./tab-ml/$run_tab
