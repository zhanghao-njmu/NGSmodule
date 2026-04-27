# Config for `ngsmodule test MethylationExtraction`.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/MethylationExtraction"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=4

Aligner=bismark
SequenceType=bisulfite
# genome_folder gets resolved from iGenomes in real runs; dry-run uses a placeholder.
