setwd("../quilt_data_for_davies/")
#dir.create("quilt_output_only_pb")
library("QUILT")

QUILT_prepare_reference(
    outputdir = "quilt_output_only_illumina",
    chr = "chr20",
    nGen = 100,
    reference_haplotype_file = "1000GP_umich_chr20_snps_genotypes_rm_dup.hap.gz",
    reference_legend_file = "1000GP_umich_chr20_snps_genotypes_rm_dup.legend.gz",
    regionStart = 60070,
    regionEnd = 64333928,
    buffer = 500000
)

system("samtools index ../quilt_data_for_davies/HG002_HiSeq30x_subsampled.filtered.sorted.illumina.2x.bwa.bam")
#system("samtools index ../quilt_data_for_davies/HG002_minimap2_hg38_sub4_sorted_AddOrReplaceReadGroups.bam")

## 
#cat(
#    "HG002_minimap2_hg38_sub4_sorted_AddOrReplaceReadGroups.bam\n",
#"HG002_HiSeq30x_subsampled.filtered.sorted.illumina.2x.bwa.bam",    
#    sep = "",
#file = "bamlist.txt")

#x <- read.table("HG002_pacbio_phasefile.txt", header = TRUE)
#x2 <- cbind(x, x)
#colnames(x2) <- c("HG002_HiSeq30x_subsampled", "HG002_pacbio_2x")
#write.table(x2, file = "phasefile.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

QUILT(
    outputdir = "quilt_output_only_illumina",
    chr = "chr20",
    regionStart = 60070,
    regionEnd = 64333928,
    buffer=500000,
    posfile = "HG002_pacbio_posfile.txt",
    phasefile = "phasefile_illumina.txt",
    bamlist = "bamlist_illumina.txt"
)
