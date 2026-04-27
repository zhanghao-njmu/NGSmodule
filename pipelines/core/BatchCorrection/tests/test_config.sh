# Config for `ngsmodule test BatchCorrection`.
# Project-scope; dry-run skips the actual sva/limma call.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/BatchCorrection"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
method=combat_seq
