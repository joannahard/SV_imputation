#!/bin/bash -l
 
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 48:00:00
#SBATCH -J sam2bam_ONT

cd /proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/ONT

module load bioinfo-tools samtools

samtools view -S -b HG002_ONT_PAD64459_Guppy_3.2_hg38.sam -@ 16 > HG002_ONT_PAD64459_Guppy_3.2.bam 2> sam2bam_ONT_log.out
