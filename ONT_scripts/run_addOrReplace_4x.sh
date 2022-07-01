#!/bin/bash -l
 
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 48:00:00
#SBATCH -J addOrReplace_ONT

module load bioinfo-tools picard

cd /proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/ONT

java -jar $PICARD_ROOT/picard.jar AddOrReplaceReadGroups \
      I=HG002_ONT_PAD64459_Guppy_3.2_sorted_4x.bam \
      O=HG002_ONT_PAD64459_Guppy_3.2_sorted_4x_AddOrReplaceReadGroups.bam \
      RGID=1 \
      RGLB=ONT \
      RGPL=ONT \
      RGPU=ONT \
      RGSM=HG002_ONT_4x
