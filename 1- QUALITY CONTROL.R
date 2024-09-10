
############################### LIBRARIES #######################
library(dada2); packageVersion("dada2") 
library(ShortRead); packageVersion("ShortRead") 
library(phyloseq); packageVersion("phyloseq") 
library(gridExtra)
require(digest)
library(readxl)
library(vegan)
library(ggplot2)
library(scales)
library(ape)
library(ggtree)
library(compositions)
library(plotly)
library(ggpubr)



#### WORKING DIRECTORY
code_path<-("C:/Users/mique/Desktop/TFM/Submarino_Giardia")
setwd(code_path)

#### Working directory to RawData 
RawDataPath <- paste(code_path,"/Giardia_Refractory_RawData_20230925_0906/",sep="")

#### Working directory to Databases
DDBBPath<-paste0(code_path,"/DDBB/")

#### Check fastq files  
list.files(RawDataPath)

# Sort files to ensures forward/reverse reads are in same order
fnFs <- sort(list.files(RawDataPath, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(RawDataPath, pattern="_R2_001.fastq"))

# Infer sample names from filenames
sample.names <- sapply(strsplit(fnFs, "_S"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(RawDataPath, fnFs)
fnRs <- file.path(RawDataPath, fnRs)

# Visualize the reads quality
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])


# Create a new file path to store filtered and trimmed reads
filt_path <- file.path(code_path, "DADA2/filtered") # place filtered files in filtered/subdirectory

# Define the name of output files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# FilterAndTrim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,verbose=T)

# How many sequences are passing the quality filters?
head(out)

plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])


# Iteratively learn an error profile from filtered reads
errF <- learnErrors(filtFs, multithread=TRUE,randomize=T,verbose=T)
errR <- learnErrors(filtRs, multithread=TRUE,randomize=T,verbose=T)

# Visualize estimated error rates 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate fastq files to speed up computation
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Inspect the dada-class object
dadaFs[[1]]
dadaRs[[1]]

# Merge denoised reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Organize ASVs into a sequence table (analogous to an OTU table)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence length
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)),main="Distribution of all PROJECT sequence length")

save.image(file = "TFM.RData")

# De novo chimera sequence detection and removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Calculate the proportion of non-chimeric ASVs
sum(seqtab.nochim)/sum(seqtab)

# Calculate number of reads obtained through each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN),rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged",  "nonchim")
rownames(track) <- sample.names
track


########################################## SEQUENCE ANALYSIS END #########################################