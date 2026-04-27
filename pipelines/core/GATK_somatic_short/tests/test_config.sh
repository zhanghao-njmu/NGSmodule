# Config for `ngsmodule test GATK_somatic_short`.
# Tumor-only mode by default.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/GATK_somatic_short"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=dna
Deduplication=TRUE
skip_bqsr=true
# normal_sample="" (unset = tumor-only mode)
