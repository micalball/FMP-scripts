
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

###################################################################################################
########################################## OBJECTS #######################################################
###################################################################################################

# ORIGINAL PHYLOSEQ
ps.clean3.def.tree
ps.melt <- psmelt(ps.clean3.def.tree)

# RELATIVE COUNTS
ps.rel <- transform_sample_counts(ps.clean3.def.tree, function(OTU) OTU/sum(OTU))
ps.melt.rel <- psmelt(ps.rel)

# GENUS
ps.genus <- tax_glom(ps.clean3.def.tree, taxrank = "Genus")
ps.melt.genus <- psmelt(ps.genus)

# GENUS RELATIVE
ps.genus.rel <- transform_sample_counts(ps.genus, function(OTU) OTU/sum(OTU))
ps.melt.genus.rel <- psmelt(ps.genus.rel)




###################################################################################################
########################################## ASVs FILTERING PHYLUM ##################################
###################################################################################################

# Aggregate by phylum
ps.rel <- transform_sample_counts(ps.clean3.def.tree, function(OTU) OTU/sum(OTU))
ps.phylum <- tax_glom(ps.clean3.def.tree, taxrank = "Phylum")
ps.phylum.rel <- transform_sample_counts(ps.phylum, function(OTU) OTU/sum(OTU))

otu_table_rel <- as.data.frame(t(otu_table(ps.phylum.rel)))
tax_table_phylum <- as.data.frame(tax_table(ps.phylum.rel))
rownames(otu_table_rel) <- tax_table_phylum$Phylum

# Step 3: Filter out taxa that have less than 0.02 relative abundance in all samples
taxa_to_keep <- rowSums(otu_table_rel >= 0.02) > 0
otu_table_filtered <- otu_table_rel[taxa_to_keep, ]

# Step 4: Also filter the tax_table to match the filtered otu_table
filtered_tax_table <- tax_table_phylum[taxa_to_keep, ]

# Step 5: Calculate the 'Others' category for each sample and add it to the OTU table
otu_table_filtered["Others", ] <- 1 - colSums(otu_table_filtered)

# Step 6: Add a row for "Others" in the filtered taxonomy table
new_tax_row <- data.frame(Kingdom = "Bacteria", 
                          Phylum = "Others", 
                          Class = NA, 
                          Order = NA, 
                          Family = NA, 
                          Genus = NA, 
                          Species = NA)
filtered_tax_table <- rbind(filtered_tax_table, new_tax_row)

# Step 7: Recreate the OTU table and Taxonomy Table objects
new_otu_table <- otu_table(as.matrix(otu_table_filtered), taxa_are_rows = TRUE)
new_tax_table <- tax_table(as.matrix(filtered_tax_table))

# Step 8: Ensure that the rownames of the OTU and taxonomy tables match
rownames(new_tax_table) <- rownames(new_otu_table)

# Step 9: Create a new phyloseq object with the filtered tables and the original sample data
ps.phylum.filtered <- phyloseq(new_otu_table, new_tax_table, sample_data(ps.phylum))

col4 <- c('royalblue', 'lightyellow',"green", 'orange', 'grey', "purple")

# Step 11: Plot the data
plot_bar(ps.phylum.filtered, x="Timepoint", fill="Phylum") + 
  facet_wrap(~Sex, scales="free_x") +
  scale_fill_manual(values = col4) +
  theme_minimal() + 
  labs(title="Relative Abundance of Taxa Grouped by Phylum",
       x="Timepoint", y="Relative Abundance")


ps.melt.phylum.rel <- psmelt(ps.phylum.filtered)
# Realizar la comparación de medias
results <- compare_means(Abundance~Condition, data = ps.melt.rel, paired = FALSE, p.adjust.method = "BH", group.by = "Phylum")
print(n= 189, results)


ggplot(ps.melt.phylum.rel, aes(x=Medication, y=Abundance*100)) +
  geom_boxplot(aes(fill=Medication),alpha=0.7) +
  geom_point() +
  scale_fill_manual(values=col) +
  theme_bw() + 
  ylab("Relative abundance") + 
  xlab("Status") +
  facet_wrap(~Phylum, scales = "free_y", nrow = 1, ncol = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=rel(1)))

col2 <- c("royalblue", "lightyellow")

ggplot(ps.melt.phylum.rel, aes(x=Condition, y=Abundance*100)) +
  geom_boxplot(aes(fill=Condition),alpha=0.7) +
  geom_point() +
  scale_fill_manual(values=col2) +
  theme_bw() + 
  ylab("Relative abundance") + 
  xlab("Status") +
  facet_wrap(~Phylum, scales = "free_y", nrow = 1, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=rel(1)))




perform_kruskal_with_posthoc_all_phyla <- function(phyloseq_object, factor_column, abundance_column) {
  # Melt the phyloseq object into a dataframe
  psmelt_df <- psmelt(phyloseq_object)
  
  # Get unique phyla
  unique_phyla <- unique(psmelt_df$Phylum)
  
  # Initialize a list to store results
  results <- list()
  
  # Iterate over each phylum
  for (phylum in unique_phyla) {
    # Subset data for the current phylum
    subset_data <- psmelt_df[psmelt_df$Phylum == phylum, ]
    
    # Perform Kruskal-Wallis test
    kruskal_result <- kruskal.test(as.formula(paste(abundance_column, "~", factor_column)), data = subset_data)
    
    # Extract the p-value
    p_value <- kruskal_result$p.value
    
    # Perform Dunn's post-hoc test
    dunn_result <- dunnTest(as.formula(paste(abundance_column, "~", factor_column)), data = subset_data, method="bh")
    
    # Store the results
    results[[phylum]] <- list(
      p_value = p_value,
      kruskal_summary = kruskal_result,
      dunn_result = dunn_result
    )
  }
  
  # Return the results for all phyla
  print(results)
  return(results)
}

# Usage example
# Apply the function to all phyla using the relative abundance dataset
all_phyla_results <- perform_kruskal_with_posthoc_all_phyla(ps.rel, "Medication", "Abundance")

# Print the results for all phyla, their p-values, and Dunn's test results
for (phylum in names(all_phyla_results)) {
  cat("Phylum:", phylum, "\n")
  cat("P-value:", all_phyla_results[[phylum]]$p_value, "\n")
  cat("Kruskal-Wallis Summary:\n")
  print(all_phyla_results[[phylum]]$kruskal_summary)
  cat("Dunn's Test Results:\n")
  print(all_phyla_results[[phylum]]$dunn_result)
  cat("\n")
}


###################################################################################################
########################################## ASVs FILTERING FAMILY ##################################
###################################################################################################

# Aggregate by family
ps.family <- tax_glom(ps.clean3.def.tree, taxrank = "Family")
ps.family.rel <- transform_sample_counts(ps.family, function(OTU) OTU/sum(OTU))

otu_table_rel <- as.data.frame(t(otu_table(ps.family.rel)))
tax_table_family <- as.data.frame(tax_table(ps.family.rel))
rownames(otu_table_rel) <- tax_table_family$Family

# Step 3: Filter out taxa that have less than 0.02 relative abundance in all samples
taxa_to_keep <- rowSums(otu_table_rel >= 0.02) > 0
otu_table_filtered <- otu_table_rel[taxa_to_keep, ]
filtered_family <- rownames(otu_table_filtered)

# Step 4: Also filter the tax_table to match the filtered otu_table
filtered_tax_table <- tax_table_family[taxa_to_keep, ]

# Step 5: Calculate the 'Others' category for each sample and add it to the OTU table
otu_table_filtered["Others", ] <- 1 - colSums(otu_table_filtered)

# Step 6: Add a row for "Others" in the filtered taxonomy table
new_tax_row <- data.frame(Kingdom = "Bacteria", 
                          Phylum = NA, 
                          Class = NA, 
                          Order = NA, 
                          Family = "Others", 
                          Genus = NA, 
                          Species = NA)
filtered_tax_table <- rbind(filtered_tax_table, new_tax_row)

# Step 7: Recreate the OTU table and Taxonomy Table objects
new_otu_table <- otu_table(as.matrix(otu_table_filtered), taxa_are_rows = TRUE)
new_tax_table <- tax_table(as.matrix(filtered_tax_table))

# Step 8: Ensure that the rownames of the OTU and taxonomy tables match
rownames(new_tax_table) <- rownames(new_otu_table)

# Step 9: Create a new phyloseq object with the filtered tables and the original sample data
ps.family.filtered <- phyloseq(new_otu_table, new_tax_table, sample_data(ps.family))

col3 <- c('tan',  'grey',  'lightyellow', "violet", 'red',"gold", 'orange','purple',"green", 'lightgreen', "yellowgreen", "black", "azure2" , 'darkred', "tomato", 'royalblue','lightpink')

# Step 11: Plot the data
plot_bar(ps.family.filtered, x="Timepoint", fill="Family") + 
  facet_wrap(~Sex, scales="free_x") +
  scale_fill_manual(values = col3) +
  theme_minimal() + 
  labs(title="Relative Abundance of Taxa Grouped by Family",
       x="Timepoint", y="Relative Abundance")

ps.melt.family.rel <- psmelt(ps.family.filtered)

# Realizar la comparación de medias
results <- compare_means(Abundance ~ Condition, data = ps.melt.rel, paired = FALSE, p.adjust.method = "BH", group.by = "Family")
print(n = 189, results)

# Boxplot for Medication
ggplot(ps.melt.family.rel, aes(x = Medication, y = Abundance * 100)) +
  geom_boxplot(aes(fill = Medication), alpha = 0.7) +
  geom_point() +
  scale_fill_manual(values = col3) +
  theme_bw() + 
  ylab("Relative abundance") + 
  xlab("Medication") +
  facet_wrap(~Family, scales = "free_y", nrow = 3, ncol = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)))

col2 <- c("royalblue", "lightyellow")

# Boxplot for Condition
ggplot(ps.melt.family.rel, aes(x = Condition, y = Abundance * 100)) +
  geom_boxplot(aes(fill = Condition), alpha = 0.7) +
  geom_point() +
  scale_fill_manual(values = col2) +
  theme_bw() + 
  ylab("Relative abundance") + 
  xlab("Condition") +
  facet_wrap(~Family, scales = "free_y", nrow = 3, ncol = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)))



perform_kruskal_with_posthoc_filtered_family <- function(phyloseq_object, factor_column, abundance_column, filtered_family) {
  # Melt the phyloseq object into a dataframe
  psmelt_df <- psmelt(phyloseq_object)
  
  # Initialize a list to store results
  all_results <- list()
  
  # Iterate over each family in the filtered list
  for (family in filtered_family) {
    # Subset data for the current family
    subset_data <- psmelt_df[psmelt_df$Family == family, ]
    
    # Perform Kruskal-Wallis test
    kruskal_result <- kruskal.test(as.formula(paste(abundance_column, "~", factor_column)), data = subset_data)
    
    # Extract the p-value
    p_value <- kruskal_result$p.value
    
    # Perform Dunn's post-hoc test
    dunn_result <- dunnTest(as.formula(paste(abundance_column, "~", factor_column)), data = subset_data, method = "bh")
    
    # Store the results
    all_results[[family]] <- list(
      p_value = p_value,
      kruskal_summary = kruskal_result,
      dunn_result = dunn_result
    )
  }
  
  # Return the results as a list
  print(all_results)
  return(all_results)
}

# Usage example
# Apply the function to the families that passed the filter using the original dataset
all_family_results <- perform_kruskal_with_posthoc_filtered_family(ps.rel, "Medication", "Abundance", filtered_family)

# Print all families, their p-values, and Dunn's test results
for (family in names(all_family_results)) {
  cat("Family:", family, "\n")
  cat("P-value:", all_family_results[[family]]$p_value, "\n")
  cat("Kruskal-Wallis Summary:\n")
  print(all_family_results[[family]]$kruskal_summary)
  cat("Dunn's Test Results:\n")
  print(all_family_results[[family]]$dunn_result)
  cat("\n")
}


###################################################################################################
########################################## ASVs FILTERING GENUS ###################################
###################################################################################################

otu_table_rel <- as.data.frame(t(otu_table(ps.genus.rel)))
tax_table_genus <- as.data.frame(tax_table(ps.genus.rel))
rownames(otu_table_rel) <- tax_table_genus$Genus

# Filter out taxa that have less than 0.02 relative abundance in all samples
taxa_to_keep <- rowSums(otu_table_rel >= 0.02) > 0
otu_table_filtered <- otu_table_rel[taxa_to_keep, ]
filtered_genera <- rownames(otu_table_filtered)

# Also filter the tax_table to match the filtered otu_table
filtered_tax_table <- tax_table_genus[taxa_to_keep, ]

# Calculate the 'Others' category for each sample and add it to the OTU table
otu_table_filtered["Others", ] <- 1 - colSums(otu_table_filtered)

# Add a row for "Others" in the filtered taxonomy table
new_tax_row <- data.frame(Kingdom = "Bacteria", 
                          Phylum = NA, 
                          Class = NA, 
                          Order = NA, 
                          Family = NA, 
                          Genus = "Others", 
                          Species = NA)
filtered_tax_table <- rbind(filtered_tax_table, new_tax_row)

# Recreate the OTU table and Taxonomy Table objects
new_otu_table <- otu_table(as.matrix(otu_table_filtered), taxa_are_rows = TRUE)
new_tax_table <- tax_table(as.matrix(filtered_tax_table))

# Ensure that the rownames of the OTU and taxonomy tables match
rownames(new_tax_table) <- rownames(new_otu_table)

# Create a new phyloseq object with the filtered tables and the original sample data
ps.filtered <- phyloseq(new_otu_table, new_tax_table, sample_data(ps.genus))

# Extract and sort the taxonomy table by Family and Genus
tax_table_genus <- as.data.frame(ps.filtered@tax_table)
tax_table_genus_sorted <- tax_table_genus[order(tax_table_genus$Family, tax_table_genus$Genus), ]

# Reorder the Genus factor in both taxonomy and OTU table
tax_table_genus_sorted$Genus <- factor(tax_table_genus_sorted$Genus, levels = unique(tax_table_genus_sorted$Genus))
otu_table_ordered <- otu_table(ps.filtered)[match(rownames(tax_table_genus_sorted), rownames(otu_table(ps.filtered))), ]

# Recreate the phyloseq object with the new ordering
ps.filtered_ordered <- phyloseq(otu_table(otu_table_ordered, taxa_are_rows = TRUE),
                                tax_table(as.matrix(tax_table_genus_sorted)),
                                sample_data(ps.filtered))

# Force the order in the plot by using `Genus` as a factor with defined levels
psmelt_genus_filt <- psmelt(ps.filtered_ordered)
psmelt_genus_filt$Genus <- factor(psmelt_genus_filt$Genus, levels = levels(tax_table_genus_sorted$Genus))

col <- c('tan',  'grey',  'lightyellow', "violet", 'red',"gold", 'orange','purple',"green", 'lightgreen', "yellowgreen", "seagreen", "azure2" , 'darkred', "tomato", 'royalblue', "darkblue", 'lightpink', "lightblue", "aquamarine","steelblue3", "salmon", "black")

# Plot ordered genus by family
ggplot(psmelt_genus_filt, aes(x=Timepoint, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity", position="stack") + 
  facet_wrap(~Sex, scales="free_x") + 
  scale_fill_manual(values = col) + 
  theme_minimal() + 
  labs(title="Relative Abundance of Taxa Grouped by Genus",
       x="Timepoint", y="Relative Abundance")

# Realizar la comparación de medias
results <- compare_means(Abundance~Condition, data = ps.melt.rel, paired = FALSE, p.adjust.method = "BH", group.by = "Genus")

# Mostrar los resultados
print(n= 189, results)

ggplot(psmelt_genus_filt, aes(x=Medication, y=Abundance*100)) +
  geom_boxplot(aes(fill=Medication),alpha=0.7) +
  geom_point() +
  scale_fill_manual(values=col3) +
  theme_bw() + 
  ylab("Relative abundance") + 
  xlab("Status") +
  facet_wrap(~Genus, scales = "free_y", nrow = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=rel(1))) 


col2 <- c("royalblue", "lightyellow")
ggplot(psmelt_genus_filt, aes(x=Condition, y=Abundance*100)) +
  geom_boxplot(aes(fill=Condition),alpha=0.7) +
  geom_point() +
  scale_fill_manual(values=col2) +
  theme_bw() + 
  ylab("Relative abundance") + 
  xlab("Status") +
  facet_wrap(~Genus, scales = "free_y", nrow = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=rel(1))) 


perform_kruskal_with_posthoc_filtered <- function(phyloseq_object, factor_column, abundance_column, genera_to_test) {
  # Melt the phyloseq object into a dataframe
  psmelt_df <- psmelt(phyloseq_object)
  
  # Initialize a list to store results
  all_results <- list()
  
  # Iterate over each genus in the filtered list
  for (genus in genera_to_test) {
    # Subset data for the current genus
    subset_data <- psmelt_df[psmelt_df$Genus == genus, ]
    
    # Perform Kruskal-Wallis test
    kruskal_result <- kruskal.test(as.formula(paste(abundance_column, "~", factor_column)), data = subset_data)
    
    # Extract the p-value
    p_value <- kruskal_result$p.value
    
    # Perform Dunn's post-hoc test
    dunn_result <- dunnTest(as.formula(paste(abundance_column, "~", factor_column)), data = subset_data, method="bh")
    
    # Store the results
    all_results[[genus]] <- list(
      p_value = p_value,
      kruskal_summary = kruskal_result,
      dunn_result = dunn_result
    )
  }
  
  # Return the results as a list
  print(all_results)
  return(all_results)
}

# Usage example
# Apply the function to the genera that passed the filter using the original dataset
all_genera_results <- perform_kruskal_with_posthoc_filtered(ps.rel, "Medication", "Abundance", filtered_genera)

# Print all genera, their p-values, and Dunn's test results
for (genus in names(all_genera_results)) {
  cat("Genus:", genus, "\n")
  cat("P-value:", all_genera_results[[genus]]$p_value, "\n")
  cat("Kruskal-Wallis Summary:\n")
  print(all_genera_results[[genus]]$kruskal_summary)
  cat("Dunn's Test Results:\n")
  print(all_genera_results[[genus]]$dunn_result)
  cat("\n")
}






