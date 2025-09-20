#!/usr/bin/env bash
set -euo pipefail
source ./config.sh
PAIRS="analysis_mutect2/00_pairs/pairs.tsv"

while IFS=$'\t' read -r patient tumor_id normal_id _; do
  [[ "$patient" == "patient" ]] && continue

  t_bam="${patient}/Tumor/03_bqsr/${tumor_id}.bqsr.bam"
  n_bam="${patient}/Normal/03_bqsr/${normal_id}.bqsr.bam"

  outdir="${patient}/Tumor/04_mutect2"
  mkdir -p "$outdir"

  args=( -R $REF_FASTA -I $t_bam -I $n_bam -tumor $tumor_id -normal $normal_id \
         --germline-resource $GERMLINE_RESOURCE -O "${outdir}/${patient}.unfiltered.vcf.gz" \
         --f1r2-tar-gz "${outdir}/${patient}.f1r2.tar.gz" \
         --native-pair-hmm-threads $THREADS )

  [[ -n "${PON_VCF:-}" ]] && args+=( --panel-of-normals $PON_VCF )

  $GATK Mutect2 "${args[@]}"
done < "$PAIRS"
