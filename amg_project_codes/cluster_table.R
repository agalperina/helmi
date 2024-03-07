library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
setwd("~/OneDrive/University of Helsinki/msc_thesis/viral_selection_2")

# Load in cluster data
helmi_clust <- read_tsv("c95_ComparisonGutCat_6.23_cluster.tsv", col_names = F)
helmi <- read_csv("subset_representatives_23.10.23.csv")
seq_6 <- read_xlsx("AVrC_6_seq.xlsx")

present <- subset(helmi_clust, X2 %in% seq_6$`catalog$contig_id`)
cluster_6_sub <- subset(helmi_clust, X1 %in% present$X1)
# keep clusters that have one of the 38 seqs and perhe
cluster_6_sub$val <- ifelse(grepl("Perhe",cluster_6_sub$X2),1,0)
cluster_6_group <- cluster_6_sub %>% 
  group_by(X1) %>% 
  summarise(val = sum(val)) %>%
  filter(val > 0) 

sequences <- cluster_6_group$X1
sequences <- gsub('Perhe1063_M_viral_212', 'GutCatV1_GPD_22771',
                  gsub('Perhe246_B4_viral_220', 'GutCatV1_GPD_100772', sequences)) # switch names from perhe ID to ID in AVrC

seq_6$HELMI <- ifelse(seq_6$`catalog$contig_id` %in% sequences,1,0) # add HELMi presence

# add taxonomy 
catalog <- read_xlsx("~/OneDrive/University of Helsinki/HUMI/phage_publication/data/unified_catalog_final.xlsx")

tax_sub <- subset(catalog,(contig_id %in% seq_6$`catalog$contig_id`)) 
tax_sub <- tax_sub %>% select(c(contig_id,genomad_Virus:phaGCN_taxonomy))
seq_6 <- cbind(seq_6,tax_sub)

# add host
host <- read_xlsx("helmi_everywhere.xlsx")
seq_6 <- cbind(seq_6,host$host)
seq_6 <- seq_6 %>% 
  rename("potential_host" = "host$host","viral_name" = "catalog$contig_id") %>%
  select(-contig_id) 

# add all samples where sequence appears
sequences <- cluster_6_group$X1
grouped_6 <- subset(helmi_clust, X1 %in% sequences)
grouped_6 <- grouped_6 %>%
  filter(grepl("Perhe",X2)) %>%
  group_by(X1) %>%
  summarise(perhe = toString(X2), .groups = 'drop') %>% 
  mutate(sample = perhe)

wanted_samples <- paste0("M|F|B4|B5|B7|B9")
grouped_6 <- grouped_6 %>% 
  mutate(sample = map(sample, ~unique(unlist(str_extract_all(.x, wanted_samples))))) 

grouped_6$X1 <- gsub('Perhe1063_M_viral_212', 'GutCatV1_GPD_22771',
                 gsub('Perhe246_B4_viral_220', 'GutCatV1_GPD_100772', grouped_6$X1))

seq_6 <- seq_6 %>% left_join(grouped_6, by = c("viral_name" = "X1"))
seq_6$sample <- gsub("NULL", NA ,seq_6$sample)

seq_6 <- seq_6[,c(1:11,21,22,12:20)]

write_xlsx(seq_6, "everywhere_table.xlsx", col_names=TRUE)
