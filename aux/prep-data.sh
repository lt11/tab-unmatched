### settings
run_tab="run-68"
pipe_type="mulo-wesone"
cd /data/ssd/lorenzot/prog/$pipe_type

### data
tar -czf train-vcf.tgz run-7/var-calls/gt/*.vcf.gz run-8/var-calls/gt/*.vcf.gz run-9/var-calls/gt/*.vcf.gz &
tar -czf valid-vcf.tgz run-10/var-calls/gt/*.vcf.gz run-12/var-calls/gt/*.vcf.gz &
tar -czf test-metmel.tgz run-11/var-calls/gt/*.vcf.gz &
tar -czf test-mixtcga.tgz run-13/var-calls/gt/*.vcf.gz run-15/var-calls/gt/*.vcf.gz &
wait
mv *tgz /data/ssd/lorenzot/prog/tab-ml

### tab-ml

cd /data/ssd/lorenzot/prog/tab-ml
mkdir -p /data/ssd/lorenzot/prog/tab-ml/$run_tab/train
mkdir -p /data/ssd/lorenzot/prog/tab-ml/$run_tab/validation
mkdir -p /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/melanoma
mkdir -p /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/mixtcga
mkdir -p /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/tnbc

cd /data/ssd/lorenzot/prog/tab-ml/
mv train-vcf.tgz /data/ssd/lorenzot/prog/tab-ml/$run_tab/train
cd /data/ssd/lorenzot/prog/tab-ml/$run_tab/train
tar -zvxf train-vcf.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf train-vcf.tgz 
gunzip *vcf.gz

cd /data/ssd/lorenzot/prog/tab-ml/
mv valid-vcf.tgz /data/ssd/lorenzot/prog/tab-ml/$run_tab/validation
cd /data/ssd/lorenzot/prog/tab-ml/$run_tab/validation
tar -zvxf valid-vcf.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf valid-vcf.tgz
gunzip *vcf.gz

cd /data/ssd/lorenzot/prog/tab-ml/
mv test-metmel.tgz /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/melanoma
cd /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/melanoma
tar -zvxf test-metmel.tgz
mv run-*/var-calls/gt/*vcf.gz .
rm -rf test-metmel.tgz
gunzip *vcf.gz

cd /data/ssd/lorenzot/prog/tab-ml/
mv test-mixtcga.tgz /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/mixtcga
cd /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/mixtcga
tar -zvxf test-mixtcga.tgz  
mv run-*/var-calls/gt/*vcf.gz .
rm -rf test-mixtcga.tgz
gunzip *vcf.gz

cd /data/ssd/lorenzot/prog/tab-ml
mv test-tnbc.tgz /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/tnbc
cd /data/ssd/lorenzot/prog/tab-ml/$run_tab/test/tnbc
tar -zvxf test-tnbc.tgz
mv tnbc/run-*/var-calls/anno-snpeff-unc/*dbs.vcf.gz .
rm -rf test-tnbc.tgz
gunzip *vcf.gz

cd /data/ssd/lorenzot/prog  
dir_tab=$(find ./tab-ml -name tab | tail -1)
cp -r ${dir_tab} ./tab-ml/$run_tab
cd ./tab-ml/$run_tab

### C
# run_tab="run-68"
# cd /Users/Lorenzo/dev/tab-ml
# scp -J eurecom -r scr/ aux/ bonette3:/data/ssd/lorenzot/prog/tab-ml/$run_tab