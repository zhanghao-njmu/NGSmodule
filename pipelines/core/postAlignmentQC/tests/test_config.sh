# Config for `ngsmodule test postAlignmentQC`.
# Auto-resolves Alignment + preAlignmentQC prerequisites (or use --no-deps).

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/postAlignmentQC"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=rna
Deduplication=NONE
# Optional: gtf=/path/to/annotation.gtf — skipped in dry-run.
