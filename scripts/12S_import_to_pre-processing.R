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
##### 9 samples (all blanks, have less than 5000 reads, this could be a cutoff)
####################################################

####################################################
# Taxonomy file
taxonomy <- read_tsv("data/raw/12S_taxonomy_BLAST96_besthit_LCA.txt") 

# Taxonomy file is missing entries for unassigned taxa, so manually refill them,
all_asvs <- colnames(ASV_table)

taxonomy_full <- tibble(`#Query` = all_asvs) %>%
  left_join(taxonomy, by = "#Query") %>%
  mutate(`#lca taxon` = replace_na(`#lca taxon`, "No data")) %>%  
  column_to_rownames("#Query") 
####################################################

####################################################
# Build phyloseq object
(UINR_12S <- phyloseq(
  otu_table(as.matrix(ASV_table), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxonomy_full)),
  sample_data(metadata)))
####################################################
# Originally, this dataset contained 1318 ASVs, 1206 of those ASVs had taxonomy information. 
# Here this visualizes the proportion of samples with taxonomy data vs no data.
taxdf <- as.data.frame(tax_table(UINR_12S))
taxdf$tax_status <- ifelse(taxdf$`#lca taxon` == "No data", "No data", "Has data")
tax_table(UINR_12S) <- tax_table(as.matrix(taxdf))
# plot
psmelt(UINR_12S) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = tax_status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("No data" = "grey70", "Has data" = "#2C7FB8")) +
  labs(x = "Sample", y = "Proportion of reads", fill = NULL) +
  theme_bw()
# Thankfully here, all samples are mostly made up of ASVs with taxonomy data.
####################################################

####################################################
# Taxonomic cleanup of the dataset
fish_taxa <- taxa_names(UINR_12S)[
  tax_table(UINR_12S)[,"#class"] %in% c("Actinopteri", "Chondrichthyes")
]
(UINR_12S_cleaned <- prune_taxa(fish_taxa, UINR_12S))
# After Taxonomy cleanup, 512 ASVs remain (of classes "Actinopteri, Chondrichthyes")
####################################################


####################################################
# Integrating OBIS data
UINR_filtered_taxonomy <- tax_table(UINR_12S_cleaned) %>%
  as("matrix") %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  as_tibble()

# OBIS taconomic occurence file
regional_OBIS <- read_csv("data/raw/OBIS_Can-Atl-Maritimes_Species_list.csv") %>%
  select(scientificName, family, genus, FBname, DemersPelag, AnaCat, Importance)

UINR_12S_with_OBIS <- UINR_filtered_taxonomy %>%
  left_join(
    regional_OBIS,
    by = c("#lca taxon" = "scientificName")
  )

write.csv(UINR_12S_with_OBIS, "data/processed/UINR_12S_w_OBIS.csv")
####################################################