#!/bin/bash

#SBATCH --array=1-16
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 15
#SBATCH -t 24:00:00
#SBATCH -J glimpse_ligate


#mkdir -p log 

module load bioinfo-tools bcftools

pwd; hostname; date;
echo "START LIGATE";

#MET=GLIMPSE_v2.0.0
#MET=GLIMPSE

PDIR=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/ligate

BIN=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/GLIMPSE_ligate_static

COV=30x_illumina_full_ref_panel #CHANGE
CHR=20

mkdir -p /proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/ligate/chr${CHR}/cov${COV}/lists

LST=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/ligate/chr${CHR}/cov${COV}/lists/list.lst
:> ${LST}

for FILE in $(ls /proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/phase/chr${CHR}/cov${COV}/imputed_*.bcf);
do
        bcftools index -f ${FILE} #Not extremely necessary, but better to do it
        echo ${FILE} >> ${LST} ##Put file names in the list
done

OUT=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/joanna_out/bcftools_illumina/ligate/chr${CHR}/cov${COV}/imputed.bcf

$BIN --input $LST --output $OUT
bcftools index -f $OUT

echo "END LIGATE";
pwd; hostname; date;
