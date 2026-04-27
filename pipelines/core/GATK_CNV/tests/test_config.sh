# Config for `ngsmodule test GATK_CNV`.
# Schema requires cnv_intervals + cnv_pon (validated up-front, even in dry-run).

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/GATK_CNV"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=dna
Deduplication=TRUE

# Schema requires these to point at real files. We use the test samples.csv
# itself as a placeholder so schema's path-syntactic check passes; the
# pipeline body short-circuits on dry-run before actually opening them.
cnv_intervals="${TEST_DIR:-.}/test_samples.csv"
cnv_pon="${TEST_DIR:-.}/test_samples.csv"

# Disable plots to avoid the .dict requirement in dry-run.
emit_plots=false
