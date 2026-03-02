# Import file
library(tidyverse)
library(phyloseq)

####################################################
# Sample metadata file
metadata <- read_csv("data/raw/UINR_metadata.csv") %>% 
  column_to_rownames("SampleName")
####################################################

####################################################
# Occurrence table file
ASV_table <- read_tsv("data/raw/12S_ASV_table.txt") %>%     # samples x taxa
  column_to_rownames("row_names")           # set sample IDs as rownames

# Inspect sample read depths
sample_sums <- rowSums(ASV_table)                          
head(sort(rowSums(ASV_table)), n = 12L) # change n to find tail end
##### 9 samples (all blanks, have less than 5000 reads, this will be cutoff)
####################################################

####################################################
# Taxonomy file
taxonomy <- read_tsv("data/raw/12S_taxonomy_BLAST96_besthit_LCA.txt") %>%  
  column_to_rownames("#Query") %>%
  filter(`#class` %in% c("Actinopteri", "Chondrichthyes", "Mammalia"))
####################################################


############################################
# Build phyloseq object

# Merging tables right away gives an error of taxa mismatch, let's investigate:
# --- Harmonize taxa (names must match) ---
taxa_keep <- intersect(colnames(ASV_table), rownames(taxonomy))
####################
# Impact tracking function:
report_taxa_differences <- function(x_names, y_names, x_label = "X", y_label = "Y") {
  only_in_x <- setdiff(x_names, y_names)
  only_in_y <- setdiff(y_names, x_names)
  
  cat(sprintf("Taxa in %s not in %s: %d\n", x_label, y_label, length(only_in_x)))
  if (length(only_in_x)) cat("  IDs:", paste(only_in_x, collapse = ", "), "\n")
  
  cat(sprintf("Taxa in %s not in %s: %d\n", y_label, x_label, length(only_in_y)))
  if (length(only_in_y)) cat("  IDs:", paste(only_in_y, collapse = ", "), "\n")
  
  invisible(list(only_in_x = only_in_x, only_in_y = only_in_y))
}

# Use it on the same intersection:
res <- report_taxa_differences(colnames(ASV_table), rownames(taxonomy),
                               x_label = "counts", y_label = "taxonomy")
###############

# Do the harmonize
ASV_table <- ASV_table[, taxa_keep, drop = FALSE]
#taxonomy   <- taxonomy[taxa_keep, , drop = FALSE]
# Then you can merge datasets :)

# Build phyloseq object
(UINR_12S <- phyloseq(
  otu_table(as.matrix(ASV_table), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metadata)))
############################################

# Originally, this dataset contained 1318 ASVs, 1206 of those ASVs had taxonomy information. 
# After Taxonomy cleanup, 558 ASVs remain (of classes "Actinopteri, Chondrichthyes, Mammalia")
# I see a good number of ASVs have less than 1000 reads, so I will apply it as a minimum threshold.  
# After read abundance filter, 397 ASVs remain. I will manually assess those against OBIS records.

UINR_filtered_taxonomy <- tax_table(UINR_12S) %>%
  as("matrix") %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  as_tibble()

# OBIS taconomic occurence file
regional_OBIS <- read_csv("data/raw/OBIS_Can-Atl-Maritimes_Species_list.csv") %>%
  select(scientificName, family, genus, FBname, DemersPelag, AnaCat, Importance)

UINR_with_OBIS <- UINR_filtered_taxonomy %>%
  left_join(
    regional_OBIS,
    by = c("#lca taxon" = "scientificName")
  )

write.csv(UINR_with_OBIS, "data/processed/UINR_w_OBIS.csv")
