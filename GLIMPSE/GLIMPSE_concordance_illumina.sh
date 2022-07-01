#!/bin/bash

#SBATCH --array=1-1
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J glimpse_concordance

echo "START CONCORDANCE";

PDIR=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/

BIN=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/GLIMPSE_concordance_static

CHR=20
COV=30x_illumina_full_ref_panel # CHANGE
#VER=2.0.0

mkdir -p /proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/concordance/chr${CHR}/cov${COV}/lists

FRQ=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/1000GP_umich_chr${CHR}_snps_sites.bcf

VAL=/proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/Illumina/data/output/HG002_HiSeq30x_subsampled/mapped/HG002_HiSeq30x_subsampled.filtered.sorted.bwa.vcf.gz

#VAL=/proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002_truth_calls/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz


EST=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/ligate/chr${CHR}/cov${COV}/imputed.bcf

LST=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/concordance/chr${CHR}/cov${COV}/lists/list.txt
:> ${LST}
echo "chr${CHR} ${FRQ} ${VAL} ${EST}" > ${LST}

OUT=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/concordance/chr${CHR}/cov${COV}/imputed

OUL=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/concordance/chr${CHR}/cov${COV}/imputed.log

$BIN --input $LST --output $OUT --log $OUL --minDP 0 --minPROB 0 --gt-validation --allele-counts --thread 1
