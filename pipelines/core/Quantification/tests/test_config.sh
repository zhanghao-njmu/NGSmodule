# Config for `ngsmodule test Quantification`.

rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/Quantification"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2
threads_featurecounts=2

Aligner=STAR
SequenceType=rna
strand_specific=0
feature_type=exon
meta_feature=gene_id
# gtf gets resolved from iGenomes in real runs; dry-run uses placeholder.
