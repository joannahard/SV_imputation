#!/bin/bash

#SBATCH --array=1-1
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J glimpse_sample

echo "START SAMPLE";

module load bioinfo-tools bcftools

EXP=bcftools_illumina
CHR=20
COV=30x_illumina_full_ref_panel
BIN=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/GLIMPSE_sample_static

mkdir -p /proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/$EXP/sample/chr${CHR}/cov${COV}

VCF=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/$EXP/ligate/chr${CHR}/cov${COV}/imputed.bcf
OUT=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/$EXP/sample/chr${CHR}/cov${COV}/phased.bcf
$BIN --input ${VCF} --solve --output ${OUT}
bcftools index -f ${OUT}
