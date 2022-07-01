#!/bin/bash -l
 
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -J sort_ONT_bam

cd /proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/ONT

module load bioinfo-tools samtools

samtools sort HG002_ONT_PAD64459_Guppy_3.2.bam -@ 16 -o HG002_ONT_PAD64459_Guppy_3.2_sorted.bam
