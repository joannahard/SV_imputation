generatePipeline <- function(infile="sample.fastq.gz", jobfile="sample.sh", samplename="sample", outdir="sample_results", reference="/proj/sens2016003/nobackup/joanna/Garrison_project/reference/GRCh38_no_alt_analysis_set.fasta", subsample=TRUE){

    samfile <- paste(outdir,"/",samplename,".sam",sep="")
    bamfileTmp <- paste(outdir,"/",samplename,"tmp.bam",sep="")
    bamfile <- paste(outdir,"/",samplename,".bam",sep="")	
    sub1bamfile <- paste(outdir,"/",samplename,"_sub1.bam",sep="")
    sub2bamfile <- paste(outdir,"/",samplename,"_sub2.bam",sep="")
    sub3bamfile <- paste(outdir,"/",samplename,"_sub3.bam",sep="")
    sub4bamfile <- paste(outdir,"/",samplename,"_sub4.bam",sep="")
    logfile1 <- paste(outdir,"/",samplename,"_log_minimap2.out",sep="")
    logfile2 <- paste(outdir,"/",samplename,"_log_samtobam.out",sep="")

    cat("#!/bin/bash -l\n", file=jobfile)
    cat(paste("#SBATCH -A sens2016003 -p core -n 16 -t 3-0 -J ",samplename,"\n",sep=""), file=jobfile, append=TRUE)

    ## Check whether all sorted bamfiles already exists

    bamsExist <- TRUE

    bamfiles <- c(bamfile)

    if(subsample == TRUE){
        bamfiles <- c(bamfile,sub1bamfile,sub2bamfile,sub3bamfile,sub4bamfile)
    }

    for(f in bamfiles){

        sortedBam <- sub(".bam","_sorted.bam",f)
        sortedBai <- sub(".bam","_sorted.bam.bai",f)

        if(FALSE %in% file.exists(c(sortedBam,sortedBai))){
            bamsExist <- FALSE
        }

    }

    if(bamsExist == FALSE){

        cat(paste("\nmodule load conda\n",sep=""), file=jobfile, append=TRUE)
        cat(paste("\nsource conda_init.sh",sep=""), file=jobfile, append=TRUE)
        cat(paste("\nconda activate /proj/sens2016003/nobackup/joanna/env/minimap2",sep=""), file=jobfile, append=TRUE)
        cat(paste("\nminimap2 -t 16 -ax map-hifi ",reference," ",infile," > ",samfile," 2> ",logfile1,"\n", sep=""), file=jobfile, append=TRUE)

        cat(paste("\nmodule load bioinfo-tools samtools picard\n",sep=""), file=jobfile, append=TRUE)

        cat(paste("\nsamtools view -S -b ",samfile," -@ 16 > ",bamfileTmp," 2> ",logfile2,"\n",sep=""), file=jobfile, append=TRUE)

        cat(paste("\nrm ",samfile,"\n",sep=""), file=jobfile, append=TRUE)

	## Add read group to bamfile

	cat(paste("\njava -jar $PICARD_ROOT/picard.jar AddOrReplaceReadGroups I=",bamfileTmp," O=",bamfile," RGID=1 RGLB=pacbio RGPL=pacbio RGPU=pacbio RGSM=",samplename,"\n",sep=""), file=jobfile, append=TRUE)	

        cat(paste("\nrm ",bamfileTmp,"\n",sep=""), file=jobfile, append=TRUE)

        ## Downsample here
        if(subsample == TRUE){

            cat(paste("\nsamtools view -@ 16 -bs 42.5 ",bamfile," > ",sub1bamfile,sep=""), file=jobfile, append=TRUE)
            cat(paste("\nsamtools view -@ 16 -bs 32.25 ",bamfile," > ",sub2bamfile,sep=""), file=jobfile, append=TRUE)
            cat(paste("\nsamtools view -@ 16 -bs 22.125 ",bamfile," > ",sub3bamfile,sep=""), file=jobfile, append=TRUE)
            cat(paste("\nsamtools view -@ 16 -bs 12.0625 ",bamfile," > ",sub4bamfile,"\n",sep=""), file=jobfile, append=TRUE)

        }

        ## Sort
        for(f in bamfiles){

            sorted <- sub(".bam","_sorted.bam",f)

            cat(paste("\nsamtools sort ",f," -@ 16 -o ",sorted,sep=""), file=jobfile, append=TRUE)

            cat(paste("\nrm ",f,"\n",sep=""), file=jobfile, append=TRUE)

        }

        cat(paste("\n",sep=""), file=jobfile, append=TRUE)

        ## Index
        for(f in bamfiles){

            sorted <- sub(".bam","_sorted.bam",f)

            cat(paste("\nsamtools index ",sorted,sep=""), file=jobfile, append=TRUE)

        }

        cat(paste("\n",sep=""), file=jobfile, append=TRUE)

    }

    ## Deepvariant

    cat(paste("\nmodule load bioinfo-tools DeepVariant\n",sep=""), file=jobfile, append=TRUE)

    cat(paste("\nMAKE_EXAMPLE_ARGS=\"variant_caller=vcf_candidate_importer,proposed_variants=/proj/sens2016003/nobackup/joanna/Garrison_project/Delaneau_collaboration/from_Delaneau/1000GP_nygc_umich_chr20_snps_only_sites_anno_maf0001.vcf.gz,vsc_min_count_snps=0,vsc_min_count_indels=0,vsc_min_fraction_snps=0,vsc_min_fraction_indels=0\"\n",sep=""), file=jobfile, append=TRUE)

    for(f in bamfiles){

        sorted <- sub(".bam","_sorted.bam",f)
        outputVCF <- sub(".bam","_deepvariant.vcf",f)
        outputgVCF <- sub(".bam","_deepvariant.gvcf",f)
        logfile <- sub(".bam","_log_deepvariant.out",f)

        ##cat(paste("\nmodule load conda\n",sep=""), file=jobfile, append=TRUE)
        ##cat(paste("\nsource conda_init.sh",sep=""), file=jobfile, append=TRUE)

        ##ulimit -u 10000
        ##module load conda
        ##source conda_init.sh
        ##module load bioinfo-tools DeepVariant

        ##mkdir -p deepvariant_forcecall_sub0
        ##mkdir -p deepvariant_forcecall_sub0/intermediate_results_dir_sub0

        cat(paste("\ndeepvariant --model_type PACBIO --ref ",reference," --reads ",sorted," --regions chr20 --output_vcf ",outputVCF," --output_gvcf ",outputgVCF," --num_shards 16 --make_examples_extra_args=\"${MAKE_EXAMPLE_ARGS}\" --postprocess_variants_extra_args=\"group_variants=false\" &> ",logfile,sep=""), file=jobfile, append=TRUE)


    }

    ## Cute SV

    cat(paste("\nmodule load bioinfo-tools conda\n",sep=""), file=jobfile, append=TRUE)
    cat(paste("\nsource conda_init.sh",sep=""), file=jobfile, append=TRUE)
    cat(paste("\nconda activate /proj/sens2016003/nobackup/joanna/env/cuteSV\n",sep=""), file=jobfile, append=TRUE)

    for(f in bamfiles){

        sorted <- sub(".bam","_sorted.bam",f)
        outputVCF <- sub(".bam","_cuteSV.vcf",f)

        cat(paste("\ncuteSV ",sorted," ",reference," ",outputVCF," . --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --threads 16 --genotype --min_support 1\n",sep=""), file=jobfile, append=TRUE)


    }

    ## BCFTools here
    cat(paste("\nmodule load bioinfo-tools bcftools\n",sep=""), file=jobfile, append=TRUE)

    vcfinput <- "/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/1000GP_umich_chr20_snps_genotypes.sites.vcf.gz"
    tsvinput <- "/proj/sens2016003/nobackup/joanna/Garrison_project/glimpse/shared/1000GP_umich_chr20_snps_genotypes.sites.tsv.gz"

    for(f in bamfiles){

        sorted <- sub(".bam","_sorted.bam",f)
        outputVCF <- sub(".bam","_bcftools_for_glimpse_full_refpanel_pacbio.vcf.gz",f)

        cat(paste("\nbcftools mpileup -f ",reference," -I -E -a 'FORMAT/DP' -T ",vcfinput," -r chr20 ",sorted," -Ou | bcftools call -Aim -C alleles -T ",tsvinput," -Oz -o ",outputVCF,"\n",sep=""), file=jobfile, append=TRUE)
	cat(paste("\nbcftools index -f ",outputVCF,"\n",sep=""), file=jobfile, append=TRUE)

   }




    ## Glimpse phase here


}

runner <- function(){

generatePipeline(infile="/proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG00438/HG00438.fastq.gz", jobfile="HG00438.sh", samplename="HG00438", outdir="HG00438_results")
generatePipeline(infile="/proj/sens2016003/nobackup/joanna/Garrison_project/data_from_wharf/HG01891/HG01891.fastq.gz", jobfile="HG01891.sh", samplename="HG01891", outdir="HG01891_results")

}




runAll <- function(infile="HPRC_samples.txt"){

samples <- read.table(infile)

   for(i in 1:nrow(samples)){

	name <- as.character(samples[i,1])
	filepath <- as.character(samples[i,2])
	
	cat(name,"\n",sep="")

	jobfile <- paste(name,".sh",sep="")
	outdir <- paste(name,"_results",sep="")

	if(!dir.exists(outdir)){
		dir.create(outdir)
	}

	generatePipeline(infile=filepath, jobfile=jobfile, samplename=name, outdir=outdir)

	cmd <- paste("sbatch ",jobfile,sep="")

	system(cmd)	

   }


}
