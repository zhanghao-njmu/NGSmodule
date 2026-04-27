# Config for `ngsmodule test GATK_CNV_PoN`.
# Project-scope; needs ≥pon_min_samples normals (Group=Normal here).

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/GATK_CNV_PoN"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
Deduplication=TRUE
pon_groups=Normal
pon_min_samples=2

# Schema requires this; use the sample CSV as a syntactic placeholder.
cnv_intervals="${TEST_DIR:-.}/test_samples.csv"
