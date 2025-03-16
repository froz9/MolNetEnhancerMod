###############################################################################
# Modified MolNetEnhancer Script with NPClassifier Integration
# -------------------------------------------------------------------------------
# This script is a modification of the original MolNetEnhancer workflow published 
# by Madeleine Ernst. In this version, NPClassifier (via the NPCTable function) 
# is used to assign chemical class annotations.
#
# The code downloads and processes GNPS network data, merges it with additional 
# annotation information, and then computes consensus chemical classes based on 
# NPClassifier results.
###############################################################################

# Load required libraries:
# - dplyr: for data manipulation
library(dplyr)
library("ChemmineR")

system('curl -d "" "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=b817262cb6114e7295fee4f73b22a3ad&view=download_cytoscape_data" -o GNPS_output_graphML.zip')
system('unzip -d GNPS_output_graphML/ GNPS_output_graphML.zip')

# Define the NAP (Network Annotation Propagation) task identifier.
nap_id = 'c4bb6b8be9e14bdebe87c6ef3abe11f6'

###############################################################################
# Determine File Paths for GNPS Files
###############################################################################
# The script checks for the existence of specific directories and files in the
# unzipped GNPS output to determine which files to load.
if ('clusterinfo_summary' %in% list.files('GNPS_output_graphML/') & 'DB_result' %in% list.files('GNPS_output_graphML/')) {
    # If standard directories exist, set file paths accordingly.
    netfile <- paste0('GNPS_output_graphML/clusterinfo_summary/',
                      list.files('GNPS_output_graphML/clusterinfo_summary/')[1])
    gnpslibfile <- paste0('GNPS_output_graphML/DB_result/',
                          list.files('GNPS_output_graphML/DB_result/')[1])
} else if ('clusterinfosummarygroup_attributes_withIDs_withcomponentID' %in% list.files('GNPS_output_graphML/')) {
    # Alternate folder naming scheme (observed in some GNPS outputs)
    netfile <- paste0('GNPS_output_graphML/clusterinfosummarygroup_attributes_withIDs_withcomponentID/',
                      list.files('GNPS_output_graphML/clusterinfosummarygroup_attributes_withIDs_withcomponentID/')[1])
    gnpslibfile <- paste0('GNPS_output_graphML/result_specnets_DB/',
                          list.files('GNPS_output_graphML/result_specnets_DB/')[1])
} else {
    # Fallback option if previous directories are not found.
    netfile <- paste0('GNPS_output_graphML/clusterinfosummary/',
                      list.files('GNPS_output_graphML/clusterinfosummary/')[1])
    gnpslibfile <- paste0('GNPS_output_graphML/result_specnets_DB/',
                          list.files('GNPS_output_graphML/result_specnets_DB/')[1])
}

###############################################################################
# Load Annotation Data from NAP and GNPS Libraries
###############################################################################
# Download NAP annotation table via URL using the defined nap_id.
nap <- read.csv(paste0("http://proteomics2.ucsd.edu/ProteoSAFe/DownloadResultFile?task=", nap_id,
                       "&block=main&file=final_out/node_attributes_table.tsv"), 
                sep = "\t", check.names = F)
# Split the SMILES by comma, take the first element, and trim whitespace
nap$ConsensusSMILES <- sapply(strsplit(nap$ConsensusSMILES, ",\\s*"), function(x) trimws(x[1]))

# Read the GNPS library file containing additional annotations.
gnpslib <- read.csv(gnpslibfile, sep = '\t', check.names = F)

# Replace various representations of empty cells with NA.
gnpslib[gnpslib == ""] <- NA
gnpslib[gnpslib == "N/A"] <- NA
gnpslib[gnpslib == " "] <- NA
gnpslib[gnpslib == "n/a"] <- NA

# Load the network file that contains the clustering information.
netfile <- read.csv(netfile, sep = '\t', check.names = F)

# Prepare Data for Merging
###############################################################################
# Select only the required columns from the network file and rename for clarity.
netfile_data <- netfile %>% 
    select(componentindex, `cluster index`) %>% 
    rename(cluster.index = `cluster index`)

# Merge GNPS library data with the network data based on cluster indices.
FinalTable_molnetenhancer <- gnpslib %>% 
    rename(cluster.index = "#Scan#") %>% 
    right_join(netfile_data, by = "cluster.index")

# Merge the NAP annotation data into the table.
FinalTable_molnetenhancer <- FinalTable_molnetenhancer %>%
    left_join(nap, by = "cluster.index")

# Create a new column 'SMILES_FINAL' that preferentially uses 'Smiles.x' unless missing,
# in which case it falls back on 'ConsensusSMILES'.
FinalTable_molnetenhancer <- FinalTable_molnetenhancer %>% 
    mutate(SMILES_FINAL = case_when(
        is.na(Smiles.x) & !is.na(ConsensusSMILES) ~ ConsensusSMILES,
        !is.na(Smiles.x) & is.na(ConsensusSMILES) ~ Smiles.x,
        !is.na(Smiles.x) & !is.na(ConsensusSMILES) ~ Smiles.x,
        TRUE ~ NA_character_
    ))
# Replace any remaining empty cells with NA.
FinalTable_molnetenhancer[FinalTable_molnetenhancer == ""] <- NA

# Filter to create a data frame with SMILES strings for NPClassifier analysis.
smiles_Classifyre <- FinalTable_molnetenhancer %>% 
    select(cluster.index, SMILES_FINAL) %>% 
    filter(!is.na(SMILES_FINAL)) %>%
    rename(smiles = SMILES_FINAL)

# smiles_Classifyre has a more rows than the required by classyfire web
# split the file in two parts
# Define batch size
BATCH_SIZE <- 999

# Split into non-overlapping chunks
split_indices <- split(1:nrow(smiles_Classifyre), ceiling(seq_along(1:nrow(smiles_Classifyre))/BATCH_SIZE))

# Write files iteratively
for (i in seq_along(split_indices)) {
  write.table(
    smiles_Classifyre[split_indices[[i]], ],
    paste0("smiles_classyfire", i, ".tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

#When the classyfire analysis finishes, download the sdf file
#https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html#Export_of_Compounds
# Create a function to process SDF files
process_sdf <- function(sdf_path) {
  sdf <- read.SDFset(sdf_path)
  blockmatrix <- datablock2ma(datablock(sdf))
  cbind(sdfid(sdf), blockmatrix) %>% 
    as.data.frame() %>% 
    setNames(c("met_ID1", colnames(blockmatrix)))
}

# Process all SDF files in a list
sdf_files <- c("12188871.sdf", "12188875.sdf")
blockmatrices <- lapply(sdf_files, process_sdf)
blockmatrix <- bind_rows(blockmatrices)


blockmatrix <- blockmatrix %>%
  mutate(across(c(Superclass, Class, Subclass), 
              ~ sub(" __ > <.*?>", "", .x)))

Classyfire_result <- blockmatrix %>% 
select(met_ID1, Superclass, 
  Class, Subclass) %>% rename ("cluster.index" = met_ID1)

rownames(Classyfire_result) <- NULL 

# Merge NPClassifier results with the original network data and rename columns to 
# indicate NPClassifier annotations.
molnetenhancer_df <- Classyfire_result %>%
  mutate(cluster.index = as.integer(cluster.index)) %>%
    right_join(netfile_data, by = "cluster.index")


# Function to retrieve the most predominant chemical class per componentindex (no changes needed here)
highestscore <- function(df, chem_level_column) {
  df %>%
    filter(
      !is.na(.data[[chem_level_column]]),
      .data[[chem_level_column]] != "",
      .data[[chem_level_column]] != "NA"
    ) %>%
    group_by(componentindex, .data[[chem_level_column]]) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    group_by(componentindex) %>%
    summarise(
      consensus_class = ifelse(n() > 0, .data[[chem_level_column]][which.max(count)], NA_character_),
      score = ifelse(n() > 0, max(count) / sum(count), NA),
      .groups = "drop"
    )
}

# Main function to assign consensus chemical classes
define_consensus_classes <- function(data) {
  # Check for required columns (no changes needed here)
  required_columns <- c("componentindex", "Superclass", "Class", "Subclass")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input data must contain the columns: componentindex, Superclass, Class, and Subclass.")
  }
  
  # Calculate consensus for each level (no changes needed here)
  superclass_consensus <- highestscore(data, "Superclass")
  class_consensus <- highestscore(data, "Class")
  subclass_consensus <- highestscore(data, "Subclass")
  
  final <- data %>%
    left_join(superclass_consensus, by = "componentindex") %>%
    rename(Superclass_Consensus = consensus_class, Superclass_Score = score) %>%
    left_join(class_consensus, by = "componentindex") %>%
    rename(Class_Consensus = consensus_class, Class_Score = score) %>%
    left_join(subclass_consensus, by = "componentindex") %>%
    rename(Subclass_Consensus = consensus_class, Subclass_Score = score) %>%
    group_by(componentindex) %>%
    mutate(
      # Propagate consensus ONLY if componentindex is NOT -1
      Superclass_Consensus = ifelse(componentindex != -1, Superclass_Consensus[1], Superclass),
      Class_Consensus = ifelse(componentindex != -1, Class_Consensus[1], Class),
      Subclass_Consensus = ifelse(componentindex != -1, Subclass_Consensus[1], Subclass)
    ) %>%
    mutate(
      # Fill ONLY if original is NA, "", or "NA" AND consensus is NOT NA
      Superclass = ifelse((is.na(Superclass) | Superclass == "" | Superclass == "NA") & !is.na(Superclass_Consensus), Superclass_Consensus, Superclass),
      Class = ifelse((is.na(Class) | Class == "" | Class == "NA") & !is.na(Class_Consensus), Class_Consensus, Class),
      Subclass = ifelse((is.na(Subclass) | Subclass == "" | Subclass == "NA") & !is.na(Subclass_Consensus), Subclass_Consensus, Subclass)
    ) %>%
    ungroup()
  
  return(final)
}

# Apply the function to define consensus classes
defined_classes <- define_consensus_classes(molnetenhancer_df) %>%
select(1,6,8,10)


# Save the final annotated table as a CSV file.
write.csv(defined_classes, file.path("Molnetenhancer.csv"), row.names = FALSE)