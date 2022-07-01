#!/bin/bash -l
 
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -J downsample_ONT

cd /proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/ONT

module load bioinfo-tools samtools


#samtools view -@ 16 -bs 42.4722 HG002_ONT_PAD64459_Guppy_3.2_sorted.bam > HG002_ONT_PAD64459_Guppy_3.2_sorted_16x.bam

#samtools view -@ 16 -bs 32.2361 HG002_ONT_PAD64459_Guppy_3.2_sorted.bam > HG002_ONT_PAD64459_Guppy_3.2_sorted_8x.bam

samtools view -@ 16 -bs 22.1180 HG002_ONT_PAD64459_Guppy_3.2_sorted.bam > HG002_ONT_PAD64459_Guppy_3.2_sorted_4x.bam

#samtools view -@ 16 -bs 12.059 HG002_ONT_PAD64459_Guppy_3.2_sorted.bam > HG002_ONT_PAD64459_Guppy_3.2_sorted_2x.bam


