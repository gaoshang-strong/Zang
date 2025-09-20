#!/usr/bin/env bash
set -euo pipefail
source ./config.sh
PAIRS="analysis_mutect2/00_pairs/pairs.tsv"

while IFS=$'\t' read -r patient tumor_id normal_id _; do
  [[ "$patient" == "patient" ]] && continue

  for sid in "$tumor_id" "$normal_id"; do
    type_label=$([[ "$sid" == "$tumor_id" ]] && echo "Tumor" || echo "Normal")
    indir="${patient}/${type_label}/01_align"
    outdir="${patient}/${type_label}/02_markdup"
    mkdir -p "$outdir"

    $GATK MarkDuplicates \
      -I "${indir}/${sid}.sorted.bam" \
      -O "${outdir}/${sid}.markdup.bam" \
      -M "${outdir}/${sid}.metrics.txt" \
      --CREATE_INDEX true
  done
done < "$PAIRS"
