# Config for `ngsmodule test Strelka2_germline`.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/Strelka2_germline"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=dna
Deduplication=TRUE
exome=FALSE
targeted=FALSE
