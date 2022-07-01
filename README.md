# SV_imputation


Low-coverage, long-read sequencing data may represent a cost-effective approach for SV discovery. To explore this possibility, we set out the develop an analytical framework for the generation of highly accurate phased haplotypes which captures the full spectrum of genetic variation, including SNPs and SVs. This is done by combining short-read whole genome sequencing data and low-coverage, long-read sequencing data.  

We begin by building haplotype scaffolds containing SNPs. In parallel, SVs are called in PacBio HiFi data derived from the same individuals. SVs are subsequently imputed into the haplotype scaffolds using MVNcall.  

We evaluate our analyses on public datasets available from the Human Pangenome Reference Consortium (HPRC), including high-coverage Illumina, PacBio HiFi and Oxford Nanopore whole genome sequencing data from 40+ individuals. To mimick low-coverage sequencing datasets, the reads in each high coverage sample are downsampled to 4 different depths of coverage, ranging from 2x to 16x.  This experimental setup allows us to evaluate imputation performance in different levels of coverage as compared to high coverage data. 


## Illumina pipeline  
This pipeline takes unaligned cram files as input and outputs processed bams and GLs in SNP sites present in a 1000GP reference panel  

Software used:  
samtools v.1.14  
bwa-mem v.0.7.17  
picard v.2.23.4  
bcftools v.1.14  

The workflow:

1. Convert cram to fastq   

   samtools collate SAMPLE.cram --reference ref.fasta -u -o SAMPLE_cram.bam SAMPLE   

   samtools fastq SAMPLE_cram.bam -1 SAMPLE_1.fq.gz.fq.gz -2 SAMPLE_2.fq.gz.fq.gz   

2. Align reads with bwa  

   bwa mem -M ref.fasta SAMPLE_1.fq.gz.fq.gz SAMPLE_2.fq.gz.fq.gz > SAMPLE.sam  

   samtools view -Sb SAMPLE.sam > SAMPLE.bam  

3. Filter reads  

   samtools view -b -q 2 -F 1028 SAMPLE.bam |  
   picard SortSam SORT_ORDER=coordinate INPUT=/dev/stdin OUTPUT=/dev/stdout |
   picard AddOrReplaceReadGroups CREATE_INDEX=true RGID=SAMPLE RGLB=SAMPLE RGPL=ILLUMINA RGSM=SAMPLE RGCN=\"NA\" RGPU=\"NA\ INPUT=/dev/stdin OUTPUT=SAMPLE_filtered.bam  
   

4. Downsample and index  

   samtools view -bs 42.5 SAMPLE_filtered.bam > SAMPLE_sub1.bam  
   samtools view -bs 32.25 SAMPLE_filtered.bam > SAMPLE_sub2.bam  
   samtools view -bs 22.125 SAMPLE_filtered.bam > SAMPLE_sub3.bam  
   samtools view -bs 12.0625 SAMPLE_filtered.bam > SAMPLE_sub4.bam  

   samtools index SAMPLE_sub1.bam  
   samtools index SAMPLE_sub2.bam  
   samtools index SAMPLE_sub3.bam  
   samtools index SAMPLE_sub4.bam   

5. Compute genotype likelihoods in reference panel sites for each bam file to be used for generating haplotype scaffolds

   bcftools mpileup -f ref_fasta -I -E -a 'FORMAT/DP' -T ref_panel.vcf -r chr20 SAMPLE_sub1.bam -Ou | bcftools call -Aim -C alleles -T ref_panel.tsv -Oz -o SAMPLE_sub1.vcf.gz  

   bcftools index -f SAMPLE_sub1.vcf.gz



## GLIMPSE pipeline (constructing haplotype scaffolds)

GLIMPSE is a phasing and imputation method for large-scale low-coverage sequencing studies which uses large reference panels.

GLIMPSE takes genotype likelihoods as input and outputs phased genotypes and genotype posteriors.


## QUILT pipeline (constructing haplotype scaffolds)

QUILT is a method for genotype imputation from low-coverage sequence using a large reference panel

QUILT takes bam files as input and outputs phased genotypes and genotype posteriors.


## PacBio pipeline  
This pipeline takes unaligned bams as input and outputs processed bams, GLs in SV calls, GLs in SNP sites present in a 1000GP reference panel


Software used:  
bam2fastq 1.3.0  
blasr v.5.3.2  
minimap2 v.2.24  
samtools v.1.14  
picard v.2.23.4  
bcftools v.1.14  
cuteSV v.1.0.13  
SURVIVOR v.1.0.7  


The workflow:  

1. Index unaligned bam and convert unaligned bam files to fastq

   pbindex SAMPLE.bam

   bam2fastq -o SAMPLE

2. Align reads with minimap2

   minimap2 -t 16 -ax map-hifi  ref.fasta SAMPLE.fastq > SAMPLE.sam

   samtools view -S -b SAMPLE.bam > SAMPLE.bam

3. Sort and index

   samtools sort SAMPLE.bam > SAMPLE_sorted.bam 

   samtools index SAMPLE_sorted.bam

4. Add or replace readgroup

   java -jar picard.jar AddOrReplaceReadGroups I=SAMPLE_sorted.bam O=SAMPLE_rg.bam RGID=1 RGLB=pacbio RGPL=pacbio RGPU=pacbio RGSM=HG002

5. Downsample and index

   samtools view -@ 16 -bs 42.5 SAMPLE_rg.bam > SAMPLE_sub1.bam
   samtools view -@ 16 -bs 32.25 SAMPLE_rg.bam > SAMPLE_sub2.bam
   samtools view -@ 16 -bs 22.125 SAMPLE_rg.bam > SAMPLE_sub3.bam
   samtools view -@ 16 -bs 12.0625 SAMPLE_rg.bam > SAMPLE_sub4.bam

   samtools index SAMPLE_sub1.bam
   samtools index SAMPLE_sub2.bam
   samtools index SAMPLE_sub3.bam
   samtools index SAMPLE_sub4.bam


6. Compute genotype likelihoods in reference panel sites for each bam file to be used for generating haplotype scaffolds

   bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ref_panel.vcf -r chr20 SAMPLE_sub1.bam -Ou | bcftools call -Aim -C alleles -T ref_panel.tsv -Oz -o SAMPLE_sub1.vcf.gz 

   bcftools index -f SAMPLE_sub1.vcf.gz


7. Call SVs in each bam file  
cuteSV SAMPLE_sub1.bam ref.fasta SAMPLE_sub1_cuteSV.vcf . --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --threads 16 --genotype --min_support 1


## Steps currently outside PacBio pipeline:  

* Perform SURVIVOR to merge every single vcf into merged.vcf  

   SURVIVOR merge cuteSV_files 1000 1 1 1 0 30 cuteSV_merged_found_in_1.vcf  


* Rerun cuteSV for each sample with -Ivcf merged.vcf (force calling step  

   cuteSV SAMPLE_sub1.bam ref.fasta SAMPLE_sub1_cuteSV_forcecall.vcf . -Ivcf cuteSV_merged_found_in_1.vcf --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --genotype --min_support 1




## ONT analysis   
This pipeline takes unaligned bams as input and outputs processed bams, GLs in SV calls, GLs in SNP sites present in a 1000GP reference panel


Software used:
bam2fastq 1.3.0
blasr v.5.3.2
.minimap2 v.2.24
.samtools v.1.14
picard v.2.23.4
bcftools v.1.14
cuteSV v.1.0.13
SURVIVOR v.1.0.7


The workflow:  

1. Align reads with minimap2  

   minimap2 -t 16 -a -z 600,200 -x map-ont SAMPLE.fastq.gz > SAMPLE.sam  

   
2. Sort and index  
   
   samtools view -S -b SAMPLE.sam > SAMPLE.bam  
  
   samtools sort SAMPLE.bam -@ 16 -o SAMPLE_sorted.bam  

3. Add or replace readgroup  

   java -jar picard.jar AddOrReplaceReadGroups I=SAMPLE_sorted.bam O=SAMPLE_rg.bam RGID=1 RGLB=ONT RGPL=ONT RGPU=ONT RGSM=HG002_ONT  


4. Downsample and index  

   samtools view -@ 16 -bs 42.5 SAMPLE_rg.bam > SAMPLE_sub1.bam  
   samtools view -@ 16 -bs 32.25 SAMPLE_rg.bam > SAMPLE_sub2.bam  
   samtools view -@ 16 -bs 22.125 SAMPLE_rg.bam > SAMPLE_sub3.bam  
   samtools view -@ 16 -bs 12.0625 SAMPLE_rg.bam > SAMPLE_sub4.bam  

   samtools index SAMPLE_sub1.bam  
   samtools index SAMPLE_sub2.bam  
   samtools index SAMPLE_sub3.bam  
   samtools index SAMPLE_sub4.bam  

5. Compute genotype likelihoods in reference panel sites for each bam file to be used for generating haplotype scaffolds  

   bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ref_panel.vcf -r chr20 SAMPLE_sub1.bam -Ou | bcftools call -Aim -C alleles -T ref_panel.tsv -Oz -o SAMPLE_sub1.vcf.gz  

   bcftools index -f SAMPLE_sub1.vcf.gz  



## Impute SVs with MVNcall

MVNcall is a program developed for genotype calling and phasing using low-coverage next-generation sequencing reads information.

MVNcall takes haplotype scaffolds as input and genotype likelihoods of polymorphic variants to be phased. MVNcall outputs a vcf file with phased polymorphic variants.

MVNcall analysis is ongoing.



For questions please email: Joanna Hård (joanna8311@gmail.com)
