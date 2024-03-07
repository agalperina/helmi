library(splitstackshape)
library(tidyverse)
library(readxl)
library(writexl)
setwd("~/OneDrive/University of Helsinki/msc_thesis/viral_selection_2")

# Determine new numbers from catalog
helmi <- read_csv("Catalogue_viralSeq_V1.03.23.csv")
quality_helmi <- helmi %>% count(checkv_quality)

#helmi <- helmi %>% separate(sample_id, into = c("Family", "Sample"), sep = "_") # use this data for plotting
#write.csv(helmi,"helmi_for_plotting.csv")

# Calculaye # OTUs in HELMi catalog
helmi_clust <- read_tsv("c95_species_vOTU_cluster.tsv", col_names = F)

# singletons - match to themselves 
helmi_group <- helmi_clust %>% group_by(X1) %>% tally()
# how many tally = 1
sum(helmi_group$n == 1) # 318018

# vOTUs present - unique seq in the left column 
length(unique(helmi_clust$X1)) # 500880
unique_OTUs <- as.data.frame(unique(helmi_clust$X1))
unique_OTUs$`unique(helmi_clust$X1)` <- gsub("viral","provirus",unique_OTUs$`unique(helmi_clust$X1)`)

#write_csv(unique_OTUs, "unique_OTU")

# how many HELMi cluster with gut catalog
helmi_cat_clust <- read_tsv("c95_ComparisonGutCat_6.23_cluster.tsv", col_names = F)
helmi_cat_clust$O2 <- ifelse(grepl("Perhe",helmi_cat_clust$X2),"helmi","catalog")

cat_clust <- helmi_cat_clust %>% group_by(X1) %>% select(-X2) %>%
  unique() %>% group_by(X1) %>% tally()

overlap <- cat_clust %>% filter(n >= 2)
sum(grepl("Perhe",overlap$X1))

avrc <- read_xlsx("unified_catalog_final.xlsx")
#gut representative seq from overlap and what other studies they are in 
avrc_gut_reps <- overlap %>% filter(!grepl("Perhe",X1))%>%
  select(X1)# just the avrc overlaps without perhe seqs

avrc_gut_seq_filter <- avrc %>% filter(avrc$contig_id %in% avrc_gut_reps$X1)%>% # keep only what is from perhe cluster
  select(GVDv1:Koonin) %>% # select the study columns
  mutate(across(everything(), ~+as.integer(.x))) %>% # change from boolean to integer
  summarise(across(where(is.numeric), sum)) # sum across columns

# change votu from perhe to gut study
avrc_perhe_reps <- helmi_cat_clust %>% 
  filter(!grepl("Perhe",X2)) %>% 
  filter(grepl("Perhe",X1)) %>% 
  select(X2)

# give us from the perhe to gut representatives in the overlap the additional catalogs those seqs cluster with 
avrc_perhe_seq_filter <- avrc %>% filter(avrc$contig_id %in% avrc_perhe_reps$X2) %>% # keep only what is from perhe cluster
  select(GVDv1:Koonin) %>% # select the study columns
  mutate(across(everything(), ~+as.integer(.x))) %>% # change from boolean to integer
  summarise(across(where(is.numeric), sum)) # sum across columns

# Combine totals
overlap_per_study <- rbind(avrc_gut_seq_filter,avrc_perhe_seq_filter)
total_overlap_per_study <- overlap_per_study %>%
  summarise(across(where(is.numeric), sum))

total_overlap_per_study <- total_overlap_per_study[,c(2,7,6,4,1,3,8,9,5)] # rearrange order
# which catalogue overlap the most with HELMi phages
overlap_prop <- rbind(total_overlap_per_study,c(20566,10021,12986,142809,33242,63424,1347,3738,54118))
overlap_prop[3,] <- (overlap_prop[2,] - overlap_prop[1,])
overlap_prop[4,] <- (overlap_prop[1,] / overlap_prop[2,])*100
rownames(overlap_prop) <- c("overlap","total","non-overlap","prop")

# 20566,10021,12986,142809,33242,63424,1347,3738,54118 - total vOTUs in each study in order of plotting_data
plotting_data <- data.frame(Study = c("CHVD","CHVD","COPSAC","COPSAC","DEVoC","DEVoC","GPD","GPD","GVDv1","GVDv1","IMGvrGut23","IMGvrGut23","Japanese4D","Japanese4D","Koonin","Koonin","MVG","MVG"),
                            Sequence = c("Overlap","Unique","Overlap","Unique","Overlap","Unique","Overlap","Unique","Overlap","Unique","Overlap","Unique","Overlap","Unique","Overlap","Unique","Overlap","Unique"),
                            Value = c(2459,18107,1945,8076,1609,11377,15201,127608,2410,30832,3364,60060,606,741,365,3373,8606,45512))

# plot
Pal2 <- c("#CA6702","#005F73")
ggplot(plotting_data, aes(fill=Sequence, y=Value, x=Study)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("Overlap with Previously Published Studies") + 
  theme_classic() + 
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), plot.title = element_text(hjust = 0.5, size = 15),
        axis.title = element_text(size = 10),legend.title = element_text(size = 10), axis.text = element_text(size = 7),
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values=Pal2)
ggsave("overlap_plot.pdf")

### HELMi in AVrC 
helmi_in_avrc <- helmi_cat_clust %>% filter(grepl("Perhe",helmi_cat_clust$X2))
test <- merge(helmi_in_avrc,helmi_clust, by.x = "X2", by.y = "X1")
test <- test %>% select(-O2,-X2.y) %>% filter(!grepl("Perhe",test$X1)) %>% unique()
