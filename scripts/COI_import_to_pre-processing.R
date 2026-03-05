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
ASV_table <- read_tsv("data/raw/COI_ASV_table.txt") %>%     # samples x taxa
  column_to_rownames("row_names")           # set sample IDs as rownames

# Inspect sample read depths
sample_sums <- rowSums(ASV_table)                          
head(sort(rowSums(ASV_table)), n = 30L) # change n to find tail end
##### 16 samples (all blanks, have less than 15000 reads, this could be a cutoff)
####################################################

####################################################
# Taxonomy file
taxonomy <- read_tsv("data/raw/COI_taxonomy_customDB-BLAST_besthit_LCA.txt") %>%  
  column_to_rownames("#Query") %>%
#  filter(`#class` %in% c("Actinopteri", "Chondrichthyes", "Mammalia"))
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
# ASVs in ASV table with no or discarded taxonomy: x ASVs

# Do the harmonize
ASV_table <- ASV_table[, taxa_keep, drop = FALSE]
#taxonomy   <- taxonomy[taxa_keep, , drop = FALSE]
# Then you can merge datasets :)

# Build phyloseq object
(UINR_COI <- phyloseq(
  otu_table(as.matrix(ASV_table), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metadata)))
############################################
