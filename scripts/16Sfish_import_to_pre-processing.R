# Import file
library(tidyverse)
library(phyloseq)
library(ggVennDiagram) # This is used at the end to compare 12S and 16S fish detections

####################################################
# Sample metadata file
metadata <- read_csv("data/raw/UINR_metadata.csv") %>% 
  column_to_rownames("SampleName")
####################################################

####################################################
# Occurrence table file
ASV_table <- read_tsv("data/raw/16S_ASV_table.txt") %>%     # samples x taxa
  column_to_rownames("row_names")           # set sample IDs as rownames

# Inspect sample read depths
sample_sums <- rowSums(ASV_table)                          
head(sort(rowSums(ASV_table)), n = 25L) # change n to find tail end
##### 18 samples (14 blanks, 4 samples have less than 10 000 reads, this could be a cutoff)
####################################################

####################################################
# Taxonomy file
taxonomy <- read_tsv("data/raw/16S_taxonomy_BLAST96_besthit_LCA.txt") %>%  
  column_to_rownames("#Query") %>%
  filter(`#class` %in% c("Actinopteri", "Chondrichthyes", "unknown class"))
####################################################
# Already quite cleaned, from X ASVs originally.
# Microbial/protist classes removed. "unkown class" contains only 1 ASV with taxonomy: Leatherback sea turtle.
# Many skates in this dataset + basking shark + spiny dogfish

############################################
# Build phyloseq object

# Merging tables right away gives an error of taxa mismatch (between ASV table and taxonomy), let's investigate:
# Ask who is in both
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
# ASVs in ASV table with no or discarded taxonomy: 239 ASVs

# --- Harmonize taxa (names must match) ---
ASV_table <- ASV_table[, taxa_keep, drop = FALSE]
#taxonomy   <- taxonomy[taxa_keep, , drop = FALSE]
# Then you can merge datasets :)

# Build phyloseq object
(UINR_16S <- phyloseq(
  otu_table(as.matrix(ASV_table), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metadata)))
############################################

# Originally, this dataset contained 632 ASVs, 401 of which had taxonomy information. 
# After taxonomy cleanup, 393 ASVs remained (of classes "Actinopteri, Chondrichthyes, unknown(seaturtle)")
# I will manually assess those against OBIS records.

UINR_filtered_taxonomy <- tax_table(UINR_16S) %>%
  as("matrix") %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  as_tibble()

# OBIS taconomic occurence file
regional_OBIS <- read_csv("data/raw/OBIS_Can-Atl-Maritimes_Species_list.csv") %>%
  select(scientificName, family, genus, FBname, DemersPelag, AnaCat, Importance)

UINR_16S_with_OBIS <- UINR_filtered_taxonomy %>%
  left_join(
    regional_OBIS,
    by = c("#lca taxon" = "scientificName")
  )

write.csv(UINR_16S_with_OBIS, "data/processed/UINR_16S_w_OBIS.csv")


UINR_16S_with_OBIS_unique <- UINR_16S_with_OBIS %>%
  distinct(`#lca taxon`, .keep_all = TRUE)
UINR_12S_with_OBIS_unique <- UINR_12S_with_OBIS %>%
  distinct(`#lca taxon`, .keep_all = TRUE)

UINR_all_fish_unique_with_OBIS <- UINR_12S_with_OBIS_unique %>%
  mutate(marker_12S = TRUE) %>%
  full_join(
    UINR_16S_with_OBIS_unique %>% mutate(marker_16S = TRUE),
    by = "#lca taxon",
    suffix = c("_12S", "_16S")
  )

UINR_all_fish_unique_with_OBIS <- UINR_all_fish_unique_with_OBIS %>%
  mutate(
    marker = case_when(
      !is.na(marker_12S) & !is.na(marker_16S) ~ "both",
      !is.na(marker_12S) ~ "12S",
      !is.na(marker_16S) ~ "16S"
    )
  )
sets <- list(
  "12S" = UINR_all_fish_unique_with_OBIS %>%
    filter(marker %in% c("12S", "both")) %>%
    pull(`#lca taxon`),

  "16S" = UINR_all_fish_unique_with_OBIS %>%
    filter(marker %in% c("16S", "both")) %>%
    pull(`#lca taxon`)
)
ggVennDiagram(sets, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#2C7FB8") +
  theme_void()


write.csv(UINR_all_fish_unique_with_OBIS, "data/processed/UINR_all_fish_unique_with_OBIS.csv")
