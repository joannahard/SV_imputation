#!/bin/bash -l
 
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 48:00:00
#SBATCH -J minimap2_hg38_ONT

cd /proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/ONT

module load conda
source conda_init.sh
conda activate /proj/sens2016003/nobackup/joanna/env/minimap2

minimap2 -t 16 -a -z 600,200 -x map-ont /proj/sens2016003/nobackup/joanna/Garrison_project/reference/GRCh38_no_alt_analysis_set.fasta HG002_ONT_PAD64459_Guppy_3.2.fastq.gz > HG002_ONT_PAD64459_Guppy_3.2_hg38.sam 2> hg38_minimap2_ONT_log.out
