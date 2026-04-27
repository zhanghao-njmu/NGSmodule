#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/igenomes.sh
#
# Path derivation for the iGenomes reference layout. Users pick a Species +
# Source + Build and the rest is computed.
#
# Standard iGenomes tree (Illumina BaseSpace convention):
#
#   $iGenomes_Dir/
#     <Species>/
#       <Source>/                    NCBI | UCSC | Ensembl
#         <Build>/                   GRCh38 | hg38 | mm10 | ...
#           Annotation/
#             Genes/genes.gtf
#             Genes/genes.bed
#             Variation/dbsnp.vcf
#           Sequence/
#             WholeGenomeFasta/genome.fa
#             BWAIndex/genome.fa
#             Bowtie2Index/genome
#             STARIndex/
#             HISAT2Index/
#             BismarkIndex/
#
# Public API:
#   ngs_resolve_igenomes <species> <source> <build>
#     Sets globals: genome, gtf, bed, dbsnp, bwa_index, bowtie2_index,
#                   star_index, hisat2_index, bismark_index, kallisto_index,
#                   salmon_index
#
# Manual override mechanism: any *_direct variable set before the call
# wins over the derived path (e.g. Genome_direct=/data/refs/custom.fa).
###############################################################################

# ---------------------------------------------------------------------------
# Helper: prefer ${name}_direct if set, else compute from $iGenomes_Dir.
#
# Args:
#   $1  variable name (e.g. genome, gtf)
#   $2  override variable name (e.g. Genome_direct, GTF_direct)
#   $3  fallback path (constructed from $iGenomes_Dir + species/source/build)
# ---------------------------------------------------------------------------
_ngs_pick() {
  local var="$1" override_name="$2" fallback="$3"
  local override="${!override_name:-}"
  if [[ -n "$override" ]]; then
    printf -v "$var" '%s' "$override"
  else
    printf -v "$var" '%s' "$fallback"
  fi
  export "${var?}"
}

# ---------------------------------------------------------------------------
# ngs_resolve_igenomes <species> <source> <build>
#
# All paths are exported so step scripts can use them directly.
# ---------------------------------------------------------------------------
ngs_resolve_igenomes() {
  local species="${1:?species required (e.g. Homo_sapiens)}"
  local source_db="${2:?source required (NCBI | UCSC | Ensembl)}"
  local build="${3:?build required (e.g. GRCh38 | hg38 | mm10)}"

  if [[ -z "${iGenomes_Dir:-}" ]]; then
    echo "ERROR: iGenomes_Dir is not set" >&2
    return 1
  fi

  local base="$iGenomes_Dir/$species/$source_db/$build"

  _ngs_pick genome           Genome_direct           "$base/Sequence/WholeGenomeFasta/genome.fa"
  _ngs_pick gtf              GTF_direct              "$base/Annotation/Genes/genes.gtf"
  _ngs_pick bed              BED_direct              "$base/Annotation/Genes/genes.bed"
  _ngs_pick dbsnp            dbSNP_direct            "$base/Annotation/Variation/dbsnp.vcf"

  _ngs_pick bwa_index        BWA_index_direct        "$base/Sequence/BWAIndex/genome.fa"
  _ngs_pick bowtie2_index    Bowtie2_index_direct    "$base/Sequence/Bowtie2Index/genome"
  _ngs_pick bowtie_index     Bowtie_index_direct     "$base/Sequence/BowtieIndex/genome"
  _ngs_pick star_index       STAR_index_direct       "$base/Sequence/STARIndex"
  _ngs_pick hisat2_index     HISAT2_index_direct     "$base/Sequence/HISAT2Index/genome"
  _ngs_pick bismark_index    Bismark_index_direct    "$base/Sequence/BismarkIndex"
  _ngs_pick kallisto_index   Kallisto_index_direct   "$base/Sequence/KallistoIndex/transcripts.idx"
  _ngs_pick salmon_index     Salmon_index_direct     "$base/Sequence/SalmonIndex"

  # Useful for downstream (e.g. variant calling):
  _ngs_pick known_snps       known_snps_direct       "$base/Annotation/Variation/known_snps.vcf.gz"
  _ngs_pick known_indels     known_indels_direct     "$base/Annotation/Variation/known_indels.vcf.gz"
}

# ---------------------------------------------------------------------------
# Resolve an aligner-specific index by name. Used by Alignment.sh which
# branches on $Aligner.
#
# Args: $1 = aligner name (bwa | bowtie2 | star | hisat2 | bismark | ...)
# Returns the index path on stdout. Empty if unknown.
# ---------------------------------------------------------------------------
ngs_aligner_index() {
  case "${1,,}" in
    bwa)              printf '%s' "${bwa_index:-}" ;;
    bwa-mem|bwa_mem)  printf '%s' "${bwa_index:-}" ;;
    bowtie)           printf '%s' "${bowtie_index:-}" ;;
    bowtie2)          printf '%s' "${bowtie2_index:-}" ;;
    star)             printf '%s' "${star_index:-}" ;;
    hisat2)           printf '%s' "${hisat2_index:-}" ;;
    bismark|bismark_bowtie2|bismark_hisat2)
                      printf '%s' "${bismark_index:-}" ;;
    kallisto)         printf '%s' "${kallisto_index:-}" ;;
    salmon)           printf '%s' "${salmon_index:-}" ;;
    *)                printf '' ;;
  esac
}

# ---------------------------------------------------------------------------
# Verify that the resolved paths exist. Prints missing items and returns
# non-zero if any are unreadable. Used at the start of pipelines to fail
# fast with a clear message.
#
# Args: list of variable names to check (e.g. ngs_check_refs genome gtf bwa_index)
# ---------------------------------------------------------------------------
ngs_check_refs() {
  local missing=0 var path
  for var in "$@"; do
    path="${!var:-}"
    if [[ -z "$path" ]]; then
      echo "ERROR: variable not set: $var" >&2
      missing=$((missing + 1))
    elif [[ ! -e "$path" ]]; then
      echo "ERROR: $var → $path (not found)" >&2
      missing=$((missing + 1))
    fi
  done
  return "$missing"
}
