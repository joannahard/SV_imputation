#!/bin/bash -l
 
#SBATCH -A sens2016003
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J bcftools_glimpse



cd /proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG002/HG002_hg38/ONT

module load bioinfo-tools bcftools

BAM=HG002_ONT_PAD64459_Guppy_3.2_sorted_2x_AddOrReplaceReadGroups.bam
VCF=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/1000GP_umich_chr20_snps_genotypes.sites.vcf.gz
TSV=/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/1000GP_umich_chr20_snps_genotypes.sites.tsv.gz
REFGEN=/proj/sens2016003/nobackup/joanna/Garrison_project/reference/GRCh38_no_alt_analysis_set.fasta
OUT=HG002_ONT_PAD64459_Guppy_3.2_sorted_2x_AddOrReplaceReadGroups_bcftools_for_glimpse_full_refpanel.vcf.gz

bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r chr20 ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}



