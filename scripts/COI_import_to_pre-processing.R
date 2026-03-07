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
##### 22 samples (all blanks, have less than 15000 reads, this could be a cutoff)
####################################################

####################################################
# Taxonomy file
taxonomy <- read_tsv("data/raw/COI_taxonomy_customDB-BLAST_besthit_LCA.txt")

# Taxonomy file is missing entries for unassigned taxa, so manually refill them,
all_asvs <- colnames(ASV_table)

taxonomy_full <- tibble(`#Query` = all_asvs) %>%
  left_join(taxonomy, by = "#Query") %>%
  mutate(`#lca taxon` = replace_na(`#lca taxon`, "No data")) %>%  
  column_to_rownames("#Query") 
####################################################

############################################
# Build phyloseq object

# --- Harmonize taxa (Make sure ASV_table and taxonomy match) ---
taxa_keep <- intersect(colnames(ASV_table), rownames(taxonomy_full))
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

# Use function on the same intersection:
res <- report_taxa_differences(colnames(ASV_table), rownames(taxonomy_full),
                               x_label = "counts", y_label = "taxonomy")
###############

# Harmonize if needed
#ASV_table <- ASV_table[, taxa_keep, drop = FALSE]
#taxonomy   <- taxonomy[taxa_keep, , drop = FALSE]
# Then you can merge datasets :)

# Build phyloseq object
(UINR_COI <- phyloseq(
  otu_table(as.matrix(ASV_table), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxonomy_full)),
  sample_data(metadata)))
############################################
# Originally, this dataset contained 15468 ASVs, 2632 of those ASVs have taxonomy information. 

# Here this visualizes the proportion of samples with taxonomy data vs no data.
taxdf <- as.data.frame(tax_table(UINR_COI))
taxdf$tax_status <- ifelse(taxdf$`#lca taxon` == "No data", "No data", "Has data")
tax_table(UINR_COI) <- tax_table(as.matrix(taxdf))
# plot
psmelt(UINR_COI) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = tax_status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("No data" = "grey70", "Has data" = "#2C7FB8")) +
  labs(x = "Sample", y = "Proportion of reads", fill = NULL) +
  theme_bw()
# It is worth noting that a large proportion of the samples are are made up of unknown sequences.
####################################################


####################################################
# Taxonomic cleanup of the dataset
tax <- as.data.frame(tax_table(UINR_COI))
### NOT WORKING NEED TO FIX





keep <- !(tax[, "#kingdom"] %in% c("Bacteria", "Pseudomonadati", "invalid taxid", "no identification") |
          tax[, "#class"]   %in% c("Actinopteri", "Chondrichthyes", "Aves") |
          tax[, "#family"]  %in% c("Bovidae", "Canidae", "Hominidae"))
# Prune phyloseq object
(UINR_COI_cleaned <- prune_taxa(keep, UINR_COI))

####################################################


# All fish detected by COI are also in the fish datasets, so they will be removed from the COI.
# Might consider removing all mammals as well, humpback whale, pacific white-sided, muskrat, castor, mole, and contaminants.


#  filter(!`#kingdom` %in% c("Bacteria", "Pseudomonadati", "invalid taxid", "no identification")) %>% 
#  filter(!`#class` %in% c("Actinopteri", "Chondrichthyes", "Aves")) %>% 
#  filter(!`#family` %in% c("Bovidae", "Canidae", "Hominidae")) 