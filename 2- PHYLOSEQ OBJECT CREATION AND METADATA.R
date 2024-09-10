
############################### LIBRARIES #######################

library(dada2); packageVersion("dada2") ## Dada2 v.1.14.0
library(ShortRead); packageVersion("ShortRead") ## ShortRead 1.44.3
library(phyloseq); packageVersion("phyloseq") ## phyloseq 1.30.0
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
library(microbiome)
library(knitr)
library("mia")
library(ecodist)    


#### WORKING DIRECTORY
code_path<-("C:/Users/mique/Desktop/TFM/Submarino_Giardia")
setwd(code_path)

############################## ASSIGN TAXA #########################################################
# TAXA ASSIGNATION
taxa<- assignTaxonomy(seqtab.nochim, paste(code_path,"/DDBB/silva_nr99_v138.1_train_set.fa.gz",sep=""), multithread=TRUE)

# SPECIES ASSIGNATION
taxa.sp <- addSpecies(taxa, paste0(code_path,"/DDBB/silva_species_assignment_v138.1.fa.gz",""),verbose=T)
save(taxa.sp, file = "taxasp.RData")

############################# TAXA CLEANING #################################

# Clean Eukaryota, Cyanobacteria and Mitochondria
taxa_clean <- taxa.sp[!taxa.sp[, 1] %in% c("Eukaryota"), ]
table(taxa_clean[, 1])
taxa_clean2 <- taxa_clean[!taxa_clean[, 2] %in% c("Cyanobacteria", "Mitochondria"), ]
table(taxa_clean2[, 2])

# Check for NA values in the PHYLUM column
any(is.na(taxa_clean2[, "Phylum"]))

# Remove rows with NA in the PHYLUM column
taxa_clean3 <- taxa_clean2[!is.na(taxa_clean2[, "Phylum"]), ]

any(is.na(taxa_clean3[, "Phylum"]))

# Check for NA values in the GENUS column
any(is.na(taxa_clean3[, "Genus"]))

# Remove rows with NA in the GENUS column
tax_table_filtered<- taxa_clean3[!is.na(taxa_clean3[, "Genus"]), ]
any(is.na(tax_table_filtered[, "Genus"]))
save(tax_table_filtered, file = "tax_table_filtered.RData")

####################################### METADATA ######################################
# Replace fasta sequences from ASV name 
taxa.silva.bkp<-tax_table_filtered
seqtab.nochim.bkp<-seqtab.nochim

# METADATA LOADING
metadata <- read_excel(paste0(code_path, "/metadata/Giardia.metadata.xls"))
metadata <- as.data.frame(metadata)

########## METADATA CONFIGURATION
# FACTORIZING VARIABLES

metadata$Sex <- factor(metadata$Sex, levels = c(1, 2),labels = c("Male", "Female"))
metadata$qPCR <- factor(metadata$qPCR, levels = c(0, 1), labels = c("Negative", "Positive"))
metadata$Condition <- factor(metadata$Condition)
metadata$Timepoint <- factor(metadata$Timepoint)
metadata$Medication_round <- factor(metadata$Medication_round)
metadata$Suplement <- factor(metadata$Suplement)
metadata$Protector <- factor(metadata$Protector)
metadata$Medication <- factor(metadata$Medication)
metadata$Sex_Condition <- factor(metadata$Sex_Condition)
metadata$Sex_Medication <- factor(metadata$Sex_Medication)
metadata$Sex_Time <- factor(metadata$Sex_Time)

# Sort samples and set row names for metadata
samples.name <- rownames(seqtab.nochim.bkp)
rownames(metadata) <- samples.name

# Check structure of metadata after setting row names
str(metadata)

# CHECK IF NAMES ARE THE SAME
sample_names_metadata <- rownames(metadata)
sample_names_seqtab <- rownames(seqtab.nochim.bkp)
all(sample_names_metadata %in% sample_names_seqtab)
all(sample_names_seqtab %in% sample_names_metadata)

# Create sample_data object
sample_data_object <- sample_data(metadata)

############################################## PHYLOSEQ CREATION ############################################

# CREATE PhyloSeq Object
ps_giardia <- phyloseq(
  otu_table(seqtab.nochim.bkp, taxa_are_rows = FALSE),
  sample_data(sample_data_object),
  tax_table(taxa.silva.bkp)
)

summary(ps_giardia)
ps_giardia
save(ps_giardia, file = "ps_giardia.RData")

# Biostrings
dna <- Biostrings::DNAStringSet(taxa_names(ps_giardia))
names(dna) <- taxa_names(ps_giardia)

# Ajuntar "ps_giardia" amb "dna"
ps <- merge_phyloseq(ps_giardia, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
any(taxa_sums(ps)==0)

ps.clean3<-prune_taxa(taxa_sums(ps)>0,ps) #remove any non-existing OTU in the phyloseq object.
any(taxa_sums(ps.clean3)==0)

# Create a data frame with taxa sums
readsumsdf <- data.frame(
  nreads = sort(taxa_sums(ps.clean3), decreasing = TRUE),
  sorted = 1:ntaxa(ps.clean3),
  type = "ASVs"
)

# Append data frame with sample sums
readsumsdf <- rbind(
  readsumsdf,
  data.frame(
    nreads = sort(sample_sums(ps.clean3), decreasing = TRUE),
    sorted = 1:nsamples(ps.clean3),
    type = "samples"
  )
)
readsumsdf$type <- factor(readsumsdf$type, levels = c("ASVs", "samples"),
                          labels = c("Reads per ASV", "Reads per sample"))





## PLOT NUMBER OF QUALITY READS
title <- "Total number of good quality reads in Giardia dataset"
ps.clean3.ReadDistribution <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) +
  geom_bar(stat = "identity") +
  ggtitle(title) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~type, nrow = 1, scales = "free")
print(ps.clean3.ReadDistribution)

## SAVE PLOT
ggsave("ps.clean3.ReadDistribution.png", plot = ps.clean3.ReadDistribution, dpi = 300)

## SAVE OBJECT
save(ps.clean3, file = "ps_clean3.RData")

############################################## CREATE PHILOGENETIC TREE ##########################################
# Export the taxonomic table to a text file
ps.clean3.tax.table <- tax_table(ps.clean3)
write.table(as.data.frame(ps.clean3.tax.table), file = "ps.clean3.tax.table2.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Export the reference sequences to a text file
ps.clean3.refseq <- refseq(ps.clean3)
write.table(as.data.frame(ps.clean3.refseq), file = "ps.clean3.refseq2.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

####### PERFOM TREE ANALYSIS 
# Read the phylogenetic tree
ps.clean3.def.refseq.tree <- read.tree("ps.clean3.def.refseq.aligned.tree")
print(ps.clean3.def.refseq.tree)

# MERGE THE PHYLOGENETIC TREE WITH PHYLOSEQ OBJECT
ps.clean3.def.tree <- merge_phyloseq(ps.clean3, ps.clean3.def.refseq.tree)

# FINAL PHYLOSEQ OBJECT
print(ps.clean3.def.tree)
summarize_phyloseq(ps.clean3.def.tree)

# Extract tree 
tree <- phy_tree(ps.clean3.def.tree)

# Plot tree
tree_plot <- ggtree(tree) +
  geom_tiplab() + # Añade etiquetas a las puntas
  theme_tree2() # Aplica un tema de árbol

ggsave("tree.png", plot = tree_plot, dpi = 300)
ggplotly(tree_plot)
save(ps.clean3.def.tree, file = "ps.clean3.def.tree.RData")



########################### RAREFACTION CURVE #############################

# check rarefaction curves ####
ps_otu_table <- data.frame(otu_table(ps.clean3.def.tree))
raremax <- min(rowSums(ps_otu_table))
col <- c('black', 'blue', 'yellow', 'red', 'orange', 'grey', 'lightpink', 'purple', 'green','darkred')
lty <- "solid"
lwd <- 2
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
out <- with(pars,
            rarecurve(ps_otu_table, step = 1000, sample = raremax, col = col,
                      lty = lty, lwd = lwd, label = FALSE))

# save plot
pdf(file.path(code_path, 'rarefaction_curve.pdf'), width = 9, height = 7)

