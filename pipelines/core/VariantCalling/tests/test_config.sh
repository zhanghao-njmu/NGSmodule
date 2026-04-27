# Config for `ngsmodule test VariantCalling`.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/VariantCalling"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=dna
Deduplication=TRUE
snp_mq=40
indel_mq=20
# genome gets resolved from iGenomes in real runs; dry-run uses placeholder.
