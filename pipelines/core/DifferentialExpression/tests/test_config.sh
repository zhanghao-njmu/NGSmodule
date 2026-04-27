# Config for `ngsmodule test DifferentialExpression`.
# Project-scope; dry-run skips the actual edgeR call.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/DifferentialExpression"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR

# Required for the contrast
group_a=control
group_b=treatment

# Defaults are baked into meta.yml params_schema; declare here for clarity.
max_padj=0.05
min_log2fc=1
min_count=10
