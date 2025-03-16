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
# - chemodiv: for cheminformatics functions, including NPCTable used for NPClassifier analysis
library(dplyr)
library(chemodiv)

###############################################################################
# Download GNPS Data
###############################################################################
# Use system calls to download the GNPS network data (in GraphML format) and unzip it.
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

###############################################################################
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
smiles_NPC <- FinalTable_molnetenhancer %>% 
  select(cluster.index, SMILES_FINAL) %>% 
  filter(!is.na(SMILES_FINAL)) %>%
  rename(smiles = SMILES_FINAL)

###############################################################################
# NPClassifier Analysis using NPCTable
###############################################################################
# Apply the NPCTable function (from the NPClassifier integration in chemodiv)
# to annotate the SMILES strings.
NPCResult <- NPCTable(smiles_NPC)  

# Select relevant columns from the NPClassifier result.
NPCResult <- NPCResult %>% select(1, 2, 3, 5, 7) 

# Merge NPClassifier results with the original network data and rename columns to 
# indicate NPClassifier annotations.
molnetenhancer_df <- NPCResult %>% 
  right_join(netfile_data, by = "cluster.index") %>% 
  rename(
    NPC_Pathway = "pathway",
    NPC_Superclass = "superclass",
    NPC_Class = "class"
  )

###############################################################################
# Function Definitions for Consensus Chemical Class Assignment
###############################################################################
# The following functions calculate the consensus chemical class for each component
# based on the NPClassifier annotations.

# Function to determine the most predominant chemical class (or level) within each component.
highestscore <- function(df, chem_level_column) {
  df %>%
    # Filter out invalid or missing chemical level entries.
    filter(
      !is.na(.data[[chem_level_column]]),
      .data[[chem_level_column]] != "",
      .data[[chem_level_column]] != "NA"
    ) %>%
    # Count the occurrences of each chemical level within each component.
    group_by(componentindex, .data[[chem_level_column]]) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    # For each component, select the chemical level with the highest count.
    group_by(componentindex) %>%
    summarise(
      consensus_class = ifelse(n() > 0, .data[[chem_level_column]][which.max(count)], NA_character_),
      score = ifelse(n() > 0, max(count) / sum(count), NA),
      .groups = "drop"
    )
}

# Main function to assign consensus chemical classes to each component.
define_consensus_classes <- function(data) {
  # Verify that all required NPClassifier columns exist in the data.
  required_columns <- c("componentindex", "NPC_Pathway", "NPC_Superclass", "NPC_Class")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input data must contain the columns: componentindex, NPC_Pathway, NPC_Superclass, and NPC_Class.")
  }
  
  # Compute consensus annotations for each chemical level.
  superclass_consensus <- highestscore(data, "NPC_Pathway")
  class_consensus <- highestscore(data, "NPC_Superclass")
  subclass_consensus <- highestscore(data, "NPC_Class")
  
  # Merge the consensus results back into the main data frame.
  final <- data %>%
    left_join(superclass_consensus, by = "componentindex") %>%
    rename(NPC_Pathway_Consensus = consensus_class, NPC_Pathway_Score = score) %>%
    left_join(class_consensus, by = "componentindex") %>%
    rename(NPC_Superclass_Consensus = consensus_class, NPC_Superclass_Score = score) %>%
    left_join(subclass_consensus, by = "componentindex") %>%
    rename(NPC_Class_Consensus = consensus_class, NPC_Class_Score = score) %>%
    group_by(componentindex) %>%
    mutate(
      # For valid components (componentindex not equal to -1), propagate the consensus annotation.
      NPC_Pathway_Consensus = ifelse(componentindex != -1, NPC_Pathway_Consensus[1], NPC_Pathway),
      NPC_Superclass_Consensus = ifelse(componentindex != -1, NPC_Superclass_Consensus[1], NPC_Superclass),
      NPC_Class_Consensus = ifelse(componentindex != -1, NPC_Class_Consensus[1], NPC_Class)
    ) %>%
    mutate(
      # Replace missing original values with the consensus if available.
      NPC_Pathway = ifelse((is.na(NPC_Pathway) | NPC_Pathway == "" | NPC_Pathway == "NA") & !is.na(NPC_Pathway_Consensus), NPC_Pathway_Consensus, NPC_Pathway),
      NPC_Superclass = ifelse((is.na(NPC_Superclass) | NPC_Superclass == "" | NPC_Superclass == "NA") & !is.na(NPC_Superclass_Consensus), NPC_Superclass_Consensus, NPC_Superclass),
      NPC_Class = ifelse((is.na(NPC_Class) | NPC_Class == "" | NPC_Class == "NA") & !is.na(NPC_Class_Consensus), NPC_Class_Consensus, NPC_Class)
    ) %>%
    ungroup()
  
  return(final)
}

###############################################################################
# Apply Consensus Function and Save Results
###############################################################################
# Calculate consensus chemical classes and select the relevant columns.
defined_classes <- define_consensus_classes(molnetenhancer_df) %>% 
  select(1, 7, 9, 11)

# Save the final annotated table as a CSV file.
write.csv(defined_classes, file.path("Molnetenhancer.csv"), row.names = FALSE)
