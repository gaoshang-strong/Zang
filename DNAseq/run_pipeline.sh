#!/usr/bin/env bash
set -euo pipefail

# Load config variables
source ./config.sh

# Create logs directory
mkdir -p analysis_mutect2/logs

echo "Starting pipeline: $(date)" >> analysis_mutect2/logs/pipeline.log

# Step 1
echo "Running 01_align_and_sort.sh ..." | tee -a analysis_mutect2/logs/pipeline.log
bash 01_align_and_sort.sh > analysis_mutect2/logs/01_align_and_sort.log 2>&1

# Step 2
echo "Running 02_markduplicates.sh ..." | tee -a analysis_mutect2/logs/pipeline.log
bash 02_markduplicates.sh > analysis_mutect2/logs/02_markduplicates.log 2>&1

# Step 3
echo "Running 03_bqsr.sh ..." | tee -a analysis_mutect2/logs/pipeline.log
bash 03_bqsr.sh > analysis_mutect2/logs/03_bqsr.log 2>&1

# Step 4
echo "Running 04_mutect2.sh ..." | tee -a analysis_mutect2/logs/pipeline.log
bash 04_mutect2.sh > analysis_mutect2/logs/04_mutect2.log 2>&1

echo "Pipeline finished: $(date)" >> analysis_mutect2/logs/pipeline.log
