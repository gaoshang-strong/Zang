#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

PAIRS="analysis_mutect2/00_pairs/pairs.tsv"

while IFS=$'\t' read -r patient tumor_id normal_id tumor_r1 tumor_r2 normal_r1 normal_r2 tumor_type normal_type; do
  [[ "$patient" == "patient" ]] && continue

  for sid in "$tumor_id" "$normal_id"; do
    type_label=$([[ "$sid" == "$tumor_id" ]] && echo "$tumor_type" || echo "$normal_type")
    outdir="${patient}/${type_label}/01_align"
    mkdir -p "$outdir"

    fq1="${FASTQ_BASE_DIR}/$( [[ "$sid" == "$tumor_id" ]] && echo "$tumor_r1" || echo "$normal_r1" )"
    fq2="${FASTQ_BASE_DIR}/$( [[ "$sid" == "$tumor_id" ]] && echo "$tumor_r2" || echo "$normal_r2" )"

    $BWA mem -t $THREADS -R "@RG\tID:${sid}\tSM:${sid}\tPL:ILLUMINA\tLB:${sid}" \
      $REF_FASTA "$fq1" "$fq2" \
      | $SAMTOOLS sort -@ $THREADS -o "${outdir}/${sid}.sorted.bam"
    $SAMTOOLS index "${outdir}/${sid}.sorted.bam"
  done
done < "$PAIRS"
