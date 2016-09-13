HOME_DIR=Watermelon/
TEST_DIR=Watermelon/tests/rules/cutadapt/end_trim_test/
TEMP_DIR=/tmp/output/
cd $TEST_DIR
cp -a $TEST_DIR/output/ $TEMP_DIR
cd $TEMP_DIR
# python can do everything above
snakemake --cores 2 \
    --snakefile $HOME_DIR/Snakefile \
    --configfile $TEST_DIR/config.yaml \
    --force 02-cutadapt/Sample_0_trimmed_R1.fastq.gz 02-cutadapt/Sample_1_trimmed_R1.fastq.gz
# python can do everything below
gunzip 02-cutadapt/*.gz
diff -r $TEST_DIR/expected/02-cutadapt 02-cutadapt