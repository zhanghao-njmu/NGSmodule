# Config for `ngsmodule test preAlignmentQC`.
#
# Designed to run cleanly in --dry-run with no real tools installed.
# For real-tool runs, drop chr22 mini-FASTQ into tests/test_data/ and
# rerun without --dry-run.

# Resolved at runtime to absolute paths by the test runner.
rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/preAlignmentQC"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=2
threads=2

# Pipeline params (validated against meta.yml params_schema).
min_quality=20
min_length=50
trim_front1=0
trim_front2=0
detect_adapter_for_pe=true
