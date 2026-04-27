# Config for `ngsmodule test CircRNA`.
# Schema requires circexplorer2_ref; use sample CSV as syntactic placeholder.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/CircRNA"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=rna
chim_segment_min=10

# Schema requires this; use the sample CSV as a syntactic placeholder
# (the pipeline body short-circuits on dry-run before actually opening it).
circexplorer2_ref="${TEST_DIR:-.}/test_samples.csv"
