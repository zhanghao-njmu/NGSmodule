# Config for `ngsmodule test GATK_germline_short`.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/GATK_germline_short"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=dna
Deduplication=TRUE

# Skip BQSR in tests since dry-run can't verify the resource bundles.
skip_bqsr=true

# Default best-practice filter expressions are pre-set in meta.yml params_schema.
# genome resolves from iGenomes in real runs; dry-run uses a placeholder.
