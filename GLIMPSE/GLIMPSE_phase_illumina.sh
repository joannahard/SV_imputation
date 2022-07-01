#!/bin/bash

#SBATCH --array=1-16
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J glimpse_phase

mkdir -p log

echo "START IMPUTATION ${MET}-${SLURM_ARRAY_TASK_ID}";
 
PDIR=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/phase

CHR=20
COV=16x_illumina_full_ref_panel # CHANGE

DATA2=$(cat /proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/chr20_chunks.txt | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1) #Read one line from the chunk file


IRG=$(echo ${DATA2} | cut -d" " -f3) #Input region
ORG=$(echo ${DATA2} | cut -d" " -f4) #Output region
printf -v IDG "%03d" ${SLURM_ARRAY_TASK_ID}

REF=1000GP_umich_chr20_snps_genotypes #Name of the reference panel
BIN=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/GLIMPSE_phase_static #Location of the binary
RP=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/${REF}.bcf #Location of the reference panel

mkdir -p /proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/phase/chr${CHR}/cov${COV} #Output directory
mkdir -p /proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/phase/chr${CHR}/cov${COV}/log #Log directory (subfolder of output)

MAP=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/chr${CHR}.b38.gmap.gz #Location of my MAP


# CHANGE
GLS=/proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/Illumina/data/output/HG002_HiSeq30x_subsampled/mapped/HG002_HiSeq30x_subsampled.filtered.sorted.illumina.16x.bwa.vcf.gz #Location of my GLS

#bam
#GLS=/proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/Illumina/data/output/HG002_HiSeq30x_subsampled/mapped/HG002_HiSeq30x_subsampled.filtered.sorted.illumina.2x.bwa.bam

OUD=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/phase/chr${CHR}/cov${COV}/imputed_${IDG}.bcf #Output name (I use IDG so it's easier to list the files in the right order later

OUL=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/phase/chr${CHR}/cov${COV}/log/impute_${IDG}.log #Location of log


#WE did all of this just to run this command: GLIMPSE_phase! 
${BIN} --input-gl ${GLS} --input-region ${IRG} --output-region ${ORG} --map ${MAP} --output ${OUD} --reference ${RP} --thread 1 --log ${OUL} 

#Run bam
#${BIN} --input-bam ${GLS} --input-region ${IRG} --output-region ${ORG} --map ${MAP} --output ${OUD} --reference ${RP} --thread 1 --log ${OUL}

echo "END IMPUTATION ${MET}-${SLURM_ARRAY_TASK_ID}";
