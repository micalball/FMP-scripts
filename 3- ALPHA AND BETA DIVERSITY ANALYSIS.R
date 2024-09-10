
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
library(microbiome)
library(knitr)
library("mia")
library(ecodist)   
library(car)  
library(FSA)    


#### WORKING DIRECTORY
code_path<-("C:/Users/mique/Desktop/TFM/Submarino_Giardia")
setwd(code_path)

###################################################################################################
############################# DIVERSITY USING ps.clean3.def.tree OBJECT ###########################
###################################################################################################

# Extract the alpha diversity measures
alpha_diversity <- estimate_richness(ps.clean3.def.tree, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
richness_data <- plot_richness(ps.clean3.def.tree, x = "Sex", measures = c("Chao1", "Shannon", "Simpson"), color = "Sex")$data
metadata$SampleID <- sample_names(ps.clean3.def.tree)
metadata <- as.data.frame(as.matrix(sample_data(ps.clean3.def.tree)))

# Extract the alpha diversity measures and merge with metadata
alpha_diversity <- estimate_richness(ps.clean3.def.tree, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
alpha_diversity <- cbind(sample_data(ps.clean3.def.tree), alpha_diversity)

######################## CHANGE THE VARIABLE TO STUDY (EXCEPT MEDICATION)
# Wilcoxon test for Chao1 
wilcox_test_chao1 <- wilcox.test(Chao1 ~ qPCR, data = alpha_diversity)
# Wilcoxon test for Observed species
wilcox_test_obs <- wilcox.test(Observed  ~ qPCR, data = alpha_diversity)
# Wilcoxon test for Shannon 
wilcox_test_shannon <- wilcox.test(Shannon ~ qPCR, data = alpha_diversity)
# Wilcoxon test for Simpson 
wilcox_test_simpson <- wilcox.test(Simpson ~ qPCR, data = alpha_diversity)

# Print results
print(wilcox_test_chao1)
print(wilcox_test_obs)
print(wilcox_test_shannon)
print(wilcox_test_simpson)


######################### MEDICATION ANALYSIS FOR ALPHA DIVERSITY
### Subset to compare the three conditions of Medication
alpha_diversity_subset <- subset(alpha_diversity, Medication %in% c("Metronidazole", "No medication"))
alpha_diversity_subset2 <- subset(alpha_diversity, Medication %in% c("Metronidazole", "Quinacrine"))
alpha_diversity_subset3 <- subset(alpha_diversity, Medication %in% c("Quinacrine", "No medication"))

# Wilcoxon test for Chao1 between Metronidazole and No medication
wilcox_test_chao1 <- wilcox.test(Chao1 ~ Medication, data = alpha_diversity_subset)
print(wilcox_test_chao1)
# Wilcoxon test for Chao1 between Metronidazole and Quinacrine
wilcox_test_chao1 <- wilcox.test(Chao1 ~ Medication, data = alpha_diversity_subset2)
print(wilcox_test_chao1)
# Wilcoxon test for Chao1 between Quinacrine and No medication
wilcox_test_chao1 <- wilcox.test(Chao1 ~ Medication, data = alpha_diversity_subset3)
print(wilcox_test_chao1)

# Wilcoxon test for Observed between Metronidazole and No medication
wilcox_test_chao1 <- wilcox.test(Observed ~ Medication, data = alpha_diversity_subset)
print(wilcox_test_chao1)
# Wilcoxon test for Observed between Metronidazole and Quinacrine
wilcox_test_chao1 <- wilcox.test(Observed ~ Medication, data = alpha_diversity_subset2)
print(wilcox_test_chao1)
# Wilcoxon test for Observed between Quinacrine and No medication
wilcox_test_chao1 <- wilcox.test(Observed ~ Medication, data = alpha_diversity_subset3)
print(wilcox_test_chao1)

# Wilcoxon test for Shannon between Metronidazole and No medication
wilcox_test_shannon <- wilcox.test(Shannon ~ Medication, data = alpha_diversity_subset)
print(wilcox_test_shannon)
# Wilcoxon test for Shannon between Metronidazole and Quinacrine
wilcox_test_shannon <- wilcox.test(Shannon ~ Medication, data = alpha_diversity_subset2)
print(wilcox_test_shannon)
# Wilcoxon test for Shannon between Quinacrine and No medication
wilcox_test_shannon <- wilcox.test(Shannon ~ Medication, data = alpha_diversity_subset3)
print(wilcox_test_shannon)

# Wilcoxon test for Simpson between Metronidazole and No medication
wilcox_test_simpson <- wilcox.test(Simpson ~ Medication, data = alpha_diversity_subset)
print(wilcox_test_simpson)
# Wilcoxon test for Simpson between Metronidazole and Quinacrine
wilcox_test_simpson <- wilcox.test(Simpson ~ Medication, data = alpha_diversity_subset2)
print(wilcox_test_simpson)
# Wilcoxon test for Simpson between Quinacrine and No medication
wilcox_test_simpson <- wilcox.test(Simpson ~ Medication, data = alpha_diversity_subset3)
print(wilcox_test_simpson)


#################################### ALPHA DIVERSITY PLOTS ####################################

alpha_diversity_df <- data.frame(SampleID = rownames(alpha_diversity), alpha_diversity, stringsAsFactors = FALSE)
metadata_df <- data.frame(SampleID = rownames(alpha_diversity), metadata$Days_post_infection, metadata$Sex, stringsAsFactors = FALSE)
combined_data <- merge(metadata_df, alpha_diversity_df, by = "SampleID")

# Reorder and take out innecessary column-
combined_data <- combined_data[order(combined_data$Sex), ]
combined_data <- combined_data[, -6]
print(combined_data)

Observed <- ggplot(combined_data, aes(x = Days_post_infection, y = Observed, color = Sex, group = Sex)) +
       geom_point(size = 3) +
       geom_line() +
       labs(title = "Observed OTUs", x = "Days Post Infection", y = "Observed OTUs") +
       theme_minimal()
ggsave("plot_observed.png", plot = Observed, width = 8, height = 6, dpi = 300)

Chao1 <- ggplot(combined_data, aes(x = Days_post_infection, y = Chao1, color = Sex, group = Sex)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "Chao1 Index", x = "Days Post Infection", y = "Chao1 Index") +
  theme_minimal()
ggsave("plot_chao1.png", plot = Chao1, width = 8, height = 6, dpi = 300)

Shannon <- ggplot(combined_data, aes(x = Days_post_infection, y = Shannon, color = Sex, group = Sex)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "Shannon Index", x = "Days Post Infection", y = "Shannon Index") +
  theme_minimal()
ggsave("plot_shannon.png", plot = Shannon, width = 8, height = 6, dpi = 300)

Simpson <- ggplot(combined_data, aes(x = Days_post_infection, y = Simpson, color = Sex, group = Sex)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "Simpson Index", x = "Days Post Infection", y = "Simpson Index") +
  theme_minimal()
ggsave("plot_simpson.png", plot = Simpson, width = 8, height = 6, dpi = 300)

##################################### Richness graphic boxplot 
# Plot richness with boxplot
ps.clean3.def.tree.diversity2 <- plot_richness(ps.clean3.def.tree, x="Sex", measures=c("Observed", "Chao1", "Shannon", "Simpson"), color="Sex") +
  geom_boxplot(aes(fill = Sex), alpha = 0.2) +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, size = 10),
        legend.position = "none")
# Display the second plot
print(ps.clean3.def.tree.diversity2)

# Save the second plot
ggsave("diversity_boxplot.png", plot = ps.clean3.def.tree.diversity2, dpi = 300)

##############################################################################################################################
############################################### NORMALITY AND HOMOSCEDASTICITY TESTS FOR PS.REL AND PS_CLR ##################################################
############################################################################################################################

# Test normality for all ASVs in ps_rel
shapiro_tests_rel <- apply(otu_table(ps.rel), 1, shapiro.test)
shapiro_pvals_rel <- sapply(shapiro_tests_rel, function(x) x$p.value)
print(shapiro_pvals_rel)

otu_table_data <- as.matrix(otu_table(ps.rel))
# Calculate the variance for each ASV across all samples
asv_variances <- apply(otu_table_data, 2, var)
print(asv_variances)

# Plot
plot(asv_variances, main = "Variance of ASVs Across All Samples", xlab = "ASVs", ylab = "Variance", pch = 19, col = "blue")



# Test normality for all ASVs in ps_clr
shapiro_tests_clr <- apply(otu_table(ps_clr), 1, shapiro.test)
shapiro_pvals_clr <- sapply(shapiro_tests_clr, function(x) x$p.value)
print(shapiro_pvals_clr)

otu_table_data <- as.matrix(otu_table(ps_clr))
# Calculate the variance for each ASV across all samples
asv_variances <- apply(otu_table_data, 2, var)  
print(asv_variances)

# Plot
plot(asv_variances, main = "Variance of ASVs Across All Samples", xlab = "ASVs", ylab = "Variance", pch = 19, col = "blue")


##########################################################################################################
########################### BETA DIVERSITY WITH WEIGHTED UNIFRAC #########################################
##########################################################################################################

# W UNIFRAC
weighted_unifrac <- phyloseq::distance(ps.clean3.def.tree, method = "wunifrac")

# PCoA using weighted UniFrac
ord_pcoa_wunifrac <- ordinate(ps.clean3.def.tree, method = "PCoA", distance = weighted_unifrac)

pcoa_complete <- plot_ordination(ps.clean3.def.tree, ord_pcoa_wunifrac, color = "Medication_round", shape = "Sex") +
  geom_point(aes(fill = Condition), size = 4, alpha = 0.5,  stroke = 1.2) + 
  scale_shape_manual(values = c(21, 22), name = "Sex", labels = c("Male", "Female")) +  
  scale_fill_manual(values = c("white", "darkgrey"), name = "Condition", labels = c("Infection", "Post-infection")) +  
  scale_color_manual(values = c("red", "orange", "purple", "green", "blue"),
                     name = "Medication", labels = c("Metronidazole R1", "Metronidazole R2", 
                                                     "No medication", "Quinacrine R1", "Quinacrine R2")) +  
  labs(
    title = "PCoA Weighted UniFrac Dissimilarity",
    x = paste("PCoA1 (", round(100 * ord_pcoa_wunifrac$values$Relative_eig[1], 2), "%)", sep = ""),
    y = paste("PCoA2 (", round(100 * ord_pcoa_wunifrac$values$Relative_eig[2], 2), "%)", sep = "")
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),  
    legend.text = element_text(size = 9)     
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21), order = 1),  
    color = guide_legend(order = 2),  
    shape = guide_legend(order = 3)   
  )
print(pcoa_complete)

result_Condition <- adonis2(weighted_unifrac ~ Condition, data = metadata, permutations = 10000)
result_Sex <- adonis2(weighted_unifrac ~ Sex, data = metadata, permutations = 10000)
result_Days <- adonis2(weighted_unifrac ~ Days_post_infection, data = metadata, permutations = 10000)
result_Timepoint <- adonis2(weighted_unifrac ~ Timepoint, data = metadata, permutations = 10000)
result_qPCR <- adonis2(weighted_unifrac ~ qPCR, data = metadata, permutations = 10000)
result_Medication_round <- adonis2(weighted_unifrac ~ Medication_round, data = metadata, permutations = 10000)
result_Medication <- adonis2(weighted_unifrac ~ Medication, data = metadata, permutations = 10000)
result_Suplement <- adonis2(weighted_unifrac ~ Suplement, data = metadata, permutations = 10000)
result_Sex_Condition <- adonis2(weighted_unifrac ~ Sex_Condition, data = metadata, permutations = 10000)
result_Sex_Medication <- adonis2(weighted_unifrac ~ Sex_Medication, data = metadata, permutations = 10000)


# p-values and R2 values from the PERMANOVA results
results_list <- list(
  Condition = result_Condition,
  Sex = result_Sex,
  Days_post_infection = result_Days,
  Timepoint = result_Timepoint,
  qPCR = result_qPCR,
  Medication_round = result_Medication_round,
  Medication = result_Medication,
  Suplement = result_Suplement,
  Sex_Condition = result_Sex_Condition,
  Sex_Medication = result_Sex_Medication
)

# Create a data frame to store p-values and R2 values
results_df <- data.frame(
  Factor = names(results_list),
  `P-value` = sapply(results_list, function(x) x$`Pr(>F)`[1]),
  `R2` = sapply(results_list, function(x) x$R2[1])
)

# Print the DataFrame
print(results_df)


#########################################################################################
################################## CLR TRANSFORMATION AND BETA DIVERSITY WITH EUCLIDEAN DISTANCES############################################
#########################################################################################

# Step 1: Convert Phyloseq object to TreeSummarizedExperiment
tse <- makeTreeSummarizedExperimentFromPhyloseq(ps.clean3.def.tree)

# Step 2: Apply CLR Transformation
tse_clr <- transformAssay(tse, assay.type = "counts", method = "clr", pseudocount = TRUE)

# Step 3: Extract the CLR-transformed data
clr_data <- assays(tse_clr)$clr
clr_data <- t(clr_data)  # Transpose to have samples in rows

# Step 4: Create a new Phyloseq object with the CLR-transformed data
ps_clr <- phyloseq(otu_table(clr_data, taxa_are_rows = FALSE), 
                   tax_table(ps.clean3.def.tree), 
                   phy_tree(ps.clean3.def.tree), 
                   sample_data(ps.clean3.def.tree))

# Step 5: Calculate Euclidean Distance on CLR-transformed data
euclidean_dist <- dist(clr_data)

# Step 6: Perform PCoA using Euclidean Distance
ord_pcoa_euclidean <- ordinate(ps_clr, method = "PCoA", distance = euclidean_dist)

# Step 7: Plot PCoA Results
pcoa_complete <- plot_ordination(ps_clr, ord_pcoa_euclidean, color = "Medication", shape = "Sex") +
  geom_point(aes(fill = Condition), size = 4, alpha = 0.5, stroke = 1.2) + 
  scale_shape_manual(values = c(21, 22), name = "Sex", labels = c("Male", "Female")) +  
  scale_fill_manual(values = c("white", "darkgrey"), name = "Condition", labels = c("Infection", "Post-infection")) +  
  scale_color_manual(values = c("tomato", "royalblue", "green"),
                     name = "Medication", labels = c("Metronidazole", "No medication", 
                                                     "Quinacrine")) +  
  labs(
    x = paste("PCoA1 (", round(100 * ord_pcoa_euclidean$values$Relative_eig[1], 2), "%)", sep = ""),
    y = paste("PCoA2 (", round(100 * ord_pcoa_euclidean$values$Relative_eig[2], 2), "%)", sep = "")
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),  
    legend.text = element_text(size = 9)     
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21), order = 1),  
    color = guide_legend(order = 2),  
    shape = guide_legend(order = 3)   
  )

# Display the plot
print(pcoa_complete)

# Step 8: Perform PERMANOVA analysis using Euclidean distances
metadata <- data.frame(sample_data(ps.clean3.def.tree))

result_Condition <- adonis2(euclidean_dist ~ Condition, data = metadata, permutations = 10000)
result_Sex <- adonis2(euclidean_dist ~ Sex, data = metadata, permutations = 10000)
result_Days <- adonis2(euclidean_dist ~ Days_post_infection, data = metadata, permutations = 10000)
result_Timepoint <- adonis2(euclidean_dist ~ Timepoint, data = metadata, permutations = 10000)
result_qPCR <- adonis2(euclidean_dist ~ qPCR, data = metadata, permutations = 10000)
result_Medication_round <- adonis2(euclidean_dist ~ Medication_round, data = metadata, permutations = 10000)
result_Medication <- adonis2(euclidean_dist ~ Medication, data = metadata, permutations = 10000)
result_Suplement <- adonis2(euclidean_dist ~ Suplement, data = metadata, permutations = 10000)
result_Sex_Condition <- adonis2(euclidean_dist ~ Sex_Condition, data = metadata, permutations = 10000)
result_Sex_Medication <- adonis2(euclidean_dist ~ Sex_Medication, data = metadata, permutations = 10000)

# Extract p-values and R2 values from the PERMANOVA results
results_list <- list(
  Condition = result_Condition,
  Sex = result_Sex,
  Days_post_infection = result_Days,
  Timepoint = result_Timepoint,
  qPCR = result_qPCR,
  Medication_round = result_Medication_round,
  Medication = result_Medication,
  Suplement = result_Suplement,
  Sex_Condition = result_Sex_Condition,
  Sex_Medication = result_Sex_Medication
)

# Create a data frame to store p-values and R2 values
results_df <- data.frame(
  Factor = names(results_list),
  `P-value` = sapply(results_list, function(x) x$`Pr(>F)`[1]),
  `R2` = sapply(results_list, function(x) x$R2[1])
)

# Print the DataFrame
print(results_df)

######################################################################################
############################# SAMPLE DISTRIBUTION WITH HEATMAP #################

beta_diversity_matrix <- as.matrix(euclidean_dist)

# Update row and column names with new labels
rownames(beta_diversity_matrix) <- metadata$Sex_Time
colnames(beta_diversity_matrix) <- metadata$Sex_Time

Heatmap <- pheatmap::pheatmap(
  beta_diversity_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  color = colorRampPalette(c("red", "green"))(50),
  main = "Euclidean Dissimilarity Heatmap"
)

ggsave("Heatmap_samples.png", plot = Heatmap, dpi = 300)
