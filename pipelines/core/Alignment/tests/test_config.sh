# Config for `ngsmodule test Alignment`.
# Auto-resolves preAlignmentQC prerequisite (or use --no-deps to skip).

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/Alignment"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=rna
sort_threads=2
