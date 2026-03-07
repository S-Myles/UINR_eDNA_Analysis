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
taxonomy <- read_tsv("data/raw/16S_taxonomy_BLAST96_besthit_LCA.txt") 

# Taxonomy file is missing entries for unassigned taxa, so manually refill them,
all_asvs <- colnames(ASV_table)

taxonomy_full <- tibble(`#Query` = all_asvs) %>%
  left_join(taxonomy, by = "#Query") %>%
  mutate(`#lca taxon` = replace_na(`#lca taxon`, "No data")) %>%  
  column_to_rownames("#Query") 
####################################################

############################################
# Build phyloseq object
(UINR_16S <- phyloseq(
  otu_table(as.matrix(ASV_table), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxonomy_full)),
  sample_data(metadata)))
############################################
# Originally, this dataset contained 632 ASVs, 401 of which had taxonomy information. 
# Here this visualizes the proportion of samples with taxonomy data vs no data.
taxdf <- as.data.frame(tax_table(UINR_16S))
taxdf$tax_status <- ifelse(taxdf$`#lca taxon` == "No data", "No data", "Has data")
tax_table(UINR_16S) <- tax_table(as.matrix(taxdf))
# plot
psmelt(UINR_16S) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = tax_status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("No data" = "grey70", "Has data" = "#2C7FB8")) +
  labs(x = "Sample", y = "Proportion of reads", fill = NULL) +
  theme_bw()
# Thankfully here, all samples are mostly made up of ASVs with taxonomy data.
####################################################

####################################################
# Taxonomic cleanup of the dataset
fish_taxa <- taxa_names(UINR_16S)[
  tax_table(UINR_16S)[,"#class"] %in% c("Actinopteri", "Chondrichthyes", "unknown class")
]
(UINR_16S_cleaned <- prune_taxa(fish_taxa, UINR_16S))
# After taxonomy cleanup, 393 ASVs remained (of classes "Actinopteri, Chondrichthyes, unknown(seaturtle)")
# Many skates in this dataset + basking shark + spiny dogfish
####################################################


####################################################
# Integrating OBIS data
UINR_filtered_taxonomy <- tax_table(UINR_16S_cleaned) %>%
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
####################################################

####################################################
# Combining all fish data from 12S and 16S assays
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
####################################################