#!/usr/bin/env bash
set -euo pipefail

# ===== 必改部分（你的环境路径）=====
# FASTQ 目录（你给的）
export FASTQ_BASE_DIR="/ShangGaoAIProjects/Zang/single_cell_data/gastric_cancer_singlecell/naturecomm2022_GSRC_multiomic/03_wes"
# 样本清单（CSV，包含你上传的列：Run title, File name 1, File name 2, type, patient）
export SAMPLE_SHEET="/ShangGaoAIProjects/Zang/single_cell_data/gastric_cancer_singlecell/naturecomm2022_GSRC_multiomic/03_wes_analysis/Sample_list.csv"

# ===== 输出根目录（遵从你的要求）=====
export WORKDIR="/ShangGaoAIProjects/Zang/single_cell_data/gastric_cancer_singlecell/naturecomm2022_GSRC_multiomic/03_wes_analysis/analysis_mutect2"
mkdir -p "$WORKDIR"/{00_pairs,01_align,02_markdup,03_bqsr,04_mutect2,05_filter,logs}

# ===== 参考与资源（按你环境替换）=====
export REF_FASTA="/ShangGaoAIProjects/tools/reference/GRCh38/BWA/genome.fa"
export GERMLINE_RESOURCE_SUB="/ShangGaoAIProjects/tools/reference/GRCh38/Mutect2_resource/af-only-gnomad.hg38.canonical.vcf.gz"
export GERMLINE_RESOURCE="/ShangGaoAIProjects/tools/reference/GRCh38/Mutect2_resource/af-only-gnomad.hg38.vcf.gz"
# Panel of Normals（可为空）
export PON_VCF="/ShangGaoAIProjects/tools/reference/GRCh38/Mutect2_resource/1000g_pon.hg38.vcf.gz"

# 置空以下变量（没有就留空）
export KNOWN_SITES_1="/ShangGaoAIProjects/tools/reference/GRCh38/Mutect2_resource/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
export KNOWN_SITES_2="/ShangGaoAIProjects/tools/reference/GRCh38/Mutect2_resource/Homo_sapiens_assembly38.dbsnp138.vcf"
export KNOWN_SITES_3="/ShangGaoAIProjects/tools/reference/GRCh38/Mutect2_resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
export TARGET_INTERVALS=""

# ===== 资源与工具 =====
export THREADS=16
export GATK="gatk"
export BWA="bwa"
export SAMTOOLS="samtools"
