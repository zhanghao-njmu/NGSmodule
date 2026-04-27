# Config for `ngsmodule test Strelka2_somatic`.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/Strelka2_somatic"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2

Aligner=STAR
SequenceType=dna
Deduplication=TRUE
exome=FALSE
targeted=FALSE
# Required for somatic calling — references the matched normal in test_samples.csv.
normal_sample=normal_sample_PE
