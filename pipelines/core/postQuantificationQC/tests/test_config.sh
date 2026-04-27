# Config for `ngsmodule test postQuantificationQC`.
# Project-scope pipeline. Dry-run skips Rscript invocation entirely.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/postQuantificationQC"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
