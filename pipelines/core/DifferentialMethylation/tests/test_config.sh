# Config for `ngsmodule test DifferentialMethylation`.
# Project-scope; dry-run skips the actual R analysis.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/DifferentialMethylation"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=bismark
group_a=control
group_b=treatment
min_cov=3
max_qval=0.05
min_meth_diff=25
