# Config for `ngsmodule test MergeCounts`.
# Project-scope pipeline: runs once across all samples in the cohort.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/MergeCounts"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
