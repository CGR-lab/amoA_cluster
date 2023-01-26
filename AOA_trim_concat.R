library(dada2)
path <- "./cutadapt_before_DADA2"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_1_001.fastq and SAMPLENAME_2_001.fastq
fnFs <- sort(list.files(path, pattern="_1_val_1.fq", full.names = TRUE)) #pattern can be modify
fnRs <- sort(list.files(path, pattern="_2_val_2.fq", full.names = TRUE)) #pattern can be modify
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not mach")
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_2_filt.fastq.gz"))
#Set truncLen and minLen according to your dataset
#AOB amoA assembly: truncLen=c(229,229), minLen = 229
#AOA amoA gap: truncLen=c(200,200), minLen = 200
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200), minLen = 200, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
out
					  