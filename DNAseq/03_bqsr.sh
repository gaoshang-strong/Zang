#!/usr/bin/env bash
set -euo pipefail
source ./config.sh
PAIRS="analysis_mutect2/00_pairs/pairs.tsv"

DO_BQSR=1
#[[ -z "${KNOWN_SITES_1:-}" || -z "${KNOWN_SITES_2:-}" ]] && DO_BQSR=0

while IFS=$'\t' read -r patient tumor_id normal_id _; do
  [[ "$patient" == "patient" ]] && continue

  for sid in "$tumor_id" "$normal_id"; do
    type_label=$([[ "$sid" == "$tumor_id" ]] && echo "Tumor" || echo "Normal")
    indir="${patient}/${type_label}/02_markdup"
    outdir="${patient}/${type_label}/03_bqsr"
    mkdir -p "$outdir"

    inbam="${indir}/${sid}.markdup.bam"
    outbam="${outdir}/${sid}.bqsr.bam"

    if [[ $DO_BQSR -eq 1 ]]; then
      table="${outdir}/${sid}.recal.table"
      $GATK BaseRecalibrator -R $REF_FASTA -I $inbam \
        --known-sites $KNOWN_SITES_1 --known-sites $KNOWN_SITES_2 --known-sites $KNOWN_SITES_3 \
        -O $table
      $GATK ApplyBQSR -R $REF_FASTA -I $inbam \
        --bqsr-recal-file $table -O $outbam --create-output-bam-index true
    else
      cp "$inbam" "$outbam"
      $SAMTOOLS index "$outbam"
    fi
  done
done < "$PAIRS"
