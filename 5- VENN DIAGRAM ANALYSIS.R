# LOAD LIBRARIES 
library(phyloseq)
library(VennDiagram)
library(grid)


# EXTRACT OTU TABLE AND TURN INTO DATAFRAME
otu_table <- as.data.frame(t(otu_table(ps.rel)))  # Transponer para que las OTUs sean columnas

# EXTRACT METADATA
metadata <- as.data.frame(as.matrix(sample_data(ps.rel)))

# CHECK
if (!all(rownames(metadata) %in% colnames(otu_table))) {
  stop("Los nombres de las muestras en los metadatos no coinciden con las columnas en la tabla de OTU.")
}


########################### VENN DIAGRAM FOR SEX_CONDITION #######################################

# LEVELS FOR SEX CONDITION
sex_levels2 <- unique(metadata$Sex_Condition)

# CREATE LIST
otu_list <- list()

# ITERATE THROUGH LEVELS
for (sex in sex_levels2) {
  samples <- rownames(metadata[metadata$Sex_Condition == sex, ])
  
  otu_subset <- otu_table[, samples, drop = FALSE]  # Asegúrate de usar drop = FALSE para mantener la matriz
  
  otus_present <- rownames(otu_subset)[rowSums(otu_subset > 0) > 0]  # Usar rowSums para identificar OTUs presentes
  
  otu_list[[sex]] <- otus_present
}

print(otu_list)

# DIAGRAM
venn.plot <- venn.diagram(
  x = otu_list,
  category.names = sex_levels2,
  filename = NULL,  # Cambiado a NULL para visualizar en lugar de guardar
  output = TRUE,
  col = "transparent",
  fill = c("red", "orange", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.20, 0.20, 0.15, 0.15)
)

# Dibujar el diagrama
grid.draw(venn.plot)


# Determine unique OTUs for each Sex_Condition category
unique_otus <- list()
for (sex_cond in names(otu_list)) {
  # Unique OTUs for this Sex_Condition
  other_conditions <- setdiff(names(otu_list), sex_cond)
  unique_otus[[sex_cond]] <- setdiff(otu_list[[sex_cond]], unlist(otu_list[other_conditions]))
}

# Extract the taxonomy table from the phyloseq object
taxonomy_table <- as.data.frame(tax_table(ps.clean3.def.tree))

# Ensure OTU IDs are correctly aligned
if (!all(rownames(taxonomy_table) %in% rownames(otu_table))) {
  stop("The OTU IDs in the taxonomy table do not match those in the OTU table.")
}

# Function to map OTUs to their Genera using the taxonomy table
get_genera <- function(otus, taxonomy) {
  unique_genera <- unique(taxonomy[otus, "Genus"])
  return(unique_genera)
}

# Create a list of unique genera for each Sex_Condition category
genera_list <- list()
for (sex_cond in names(unique_otus)) {
  genera_list[[sex_cond]] <- get_genera(unique_otus[[sex_cond]], taxonomy_table)
}

# Print the unique genera for each Sex_Condition category
for (sex_cond in names(genera_list)) {
  cat("\nUnique Genera for", sex_cond, ":\n")
  print(genera_list[[sex_cond]])
}

############################################# COMMON ASV COMPARISION 

# Find ASVs common to all groups
common_all_groups <- Reduce(intersect, otu_list)

# Map these ASVs to their corresponding genera
genera_common_all_groups <- get_genera(common_all_groups, taxonomy_table)

# Print the results
cat("\nASVs common to all groups:\n")
print(common_all_groups)

cat("\nGenera common to all groups:\n")
print(genera_common_all_groups)

############################################# POST INFECTION COMPARISION 

# Find ASVs common to post-infection groups
common_post_infection <- intersect(otu_list[["Male post infection"]], otu_list[["Female post infection"]])

# Exclude ASVs that are also in the infected groups
common_post_infection_only <- setdiff(common_post_infection, unlist(otu_list[c("Male infected", "Female infected")]))

# Map these ASVs to their corresponding genera
genera_common_post_infection_only <- get_genera(common_post_infection_only, taxonomy_table)

# Print the results
cat("\nASVs common to post-infection only:\n")
print(common_post_infection_only)

cat("\nGenera common to post-infection only:\n")
print(genera_common_post_infection_only)

############################################# INFECTION COMPARISION 

# Find ASVs common to infection groups
common_infection <- intersect(otu_list[["Male infected"]], otu_list[["Female infected"]])

# Exclude ASVs that are also in the post-infection groups
common_infection_only <- setdiff(common_infection, unlist(otu_list[c("Male post infection", "Female post infection")]))

# Map these ASVs to their corresponding genera
genera_common_infection_only <- get_genera(common_infection_only, taxonomy_table)

# Print the results
cat("\nASVs common to infection only:\n")
print(common_infection_only)

cat("\nGenera common to infection only:\n")
print(genera_common_infection_only)



########################################################################################################
######################################### VENN DIAGRAM FOR MEDICATION ##################################
########################################################################################################

# MEDICATION LEVELS
medication_levels <- unique(metadata$Medication)

# PREPARE LIST
otu_list <- list()

#ITERATE THROUGH MEDICATION LEVELS
for (med in medication_levels) {
  samples <- rownames(metadata[metadata$Medication == med, ])
  
  otu_subset <- otu_table[, samples, drop = FALSE]  
  
  otus_present <- rownames(otu_subset)[rowSums(otu_subset > 0) > 0]  
  
  otu_list[[med]] <- otus_present
}

# ASV LIST
print(otu_list)

#VENN DIAGRAM 
venn.plot <- venn.diagram(
  x = otu_list,
  category.names = medication_levels,
  filename = NULL,  
  output = TRUE,
  col = "transparent",
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.06, 0.06, 0.09)
)

# PLOT
grid.draw(venn.plot)

# ASV FOR EACH MEDICATION CATHEGORY
unique_otus <- list()
for (med in names(otu_list)) {
  # UNIQUE OTUS FOR SPECIFIC MEDICATION
  other_meds <- setdiff(names(otu_list), med)
  unique_otus[[med]] <- setdiff(otu_list[[med]], unlist(otu_list[other_meds]))
}

# EXTRACT TAXONOMY TABLE
taxonomy_table <- as.data.frame(tax_table(ps.rel))

# CHECK
if (!all(rownames(taxonomy_table) %in% rownames(otu_table))) {
  stop("THE IDs OF ASV IN THE TAXONOMY TABLE DO NOT COINCIDE WITH ASV TABLE.")
}

# MAP ASV TO EACH GENERA BASED ON TAXONOMY TABLE
get_genera <- function(otus, taxonomy) {
  unique_genera <- unique(taxonomy[otus, "Genus"])
  return(unique_genera)
}

# LIST OF UNIQUE GENERA FOR EACH MEDICATION
genera_list <- list()
for (med in names(unique_otus)) {
  genera_list[[med]] <- get_genera(unique_otus[[med]], taxonomy_table)
}

# PRINT
for (med in names(genera_list)) {
  cat("\nUnique Genera for", med, ":\n")
  print(genera_list[[med]])
        print(unique_otus[[med]])
}



############################################# COMMON TO ALL COMPARISION

# Find ASVs common to all groups
common_all_groups <- Reduce(intersect, otu_list)

# Map these ASVs to their corresponding genera
genera_common_all_groups <- get_genera(common_all_groups, taxonomy_table)

# Print the results
cat("\nASVs common to all groups:\n")
print(common_all_groups)

cat("\nGenera common to all groups:\n")
print(genera_common_all_groups)

############################################# MEDICATED COMPARISION

# Assuming 'otu_list' contains lists of ASVs for each group (e.g., "Metronidazole", "Quinacrine", and "No medication")
# Find ASVs common to Metronidazole and Quinacrine post-infection groups
common_medication <- intersect(otu_list[["Metronidazole"]], otu_list[["Quinacrine"]])

# Exclude ASVs that are also present in the "No medication" group
common_medication_only <- setdiff(common_medication, unlist(otu_list[["No medication"]]))

# Map these ASVs to their corresponding genera using the taxonomy table
genera_common_medication_only <- get_genera(common_medication_only, taxonomy_table)

# Print the results
cat("\nASVs common to post-infection groups (excluding No medication group):\n")
print(common_medication_only)

cat("\nGenera common to post-infection groups (excluding No medication group):\n")
print(genera_common_medication_only)
