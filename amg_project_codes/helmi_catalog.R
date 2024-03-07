library(splitstackshape)
library(readxl)
library(writexl)
library(gridExtra)
library(tidyverse)

setwd("~/OneDrive/University of Helsinki/msc_thesis/viral_selection_2")

helmi <- read_csv("subset_representatives_23.10.23.csv")
phatyp <- read_csv("Phatyp_HelmiVirus.tsv")
phagcn <- read_csv("PhaGCN_HelmiVirus.tsv")
genomad <- read_tsv("Genomad_HelmiVirus.tsv")

# Merge genomad
genomad <- genomad %>%
  select(c("seq_name","topology","taxonomy")) %>%
  rename("viral_name" = "seq_name","genomad_taxonomy" = "taxonomy")
genomad$viral_name <- gsub("\\|.*","",genomad$viral_name)

genomad_sub <- subset(genomad, viral_name %in% helmi$viral_name)
genomad_missing <- subset(helmi, !(viral_name %in% genomad$viral_name)) # not present in genomad but present in helmi
helmi <- merge(helmi, genomad_sub, all=T)
# remove duplicates
helmi <- helmi[!duplicated(helmi$viral_name),] 

# Remove plasmids - NA in genomad merge, split taxonomy
# plasmids
plasmids <- helmi %>% filter(is.na(topology))

helmi <- helmi %>% 
  select(-c("sample_id","contig_id","gene_count":"host_genes","completeness":"kmer_freq","Avg_fold":"Covered_bases")) %>%
  filter(!is.na(topology)) %>%
  separate(genomad_taxonomy, c("genomad_Virus", "genomad_Realm", "genomad_Kingdom","genomad_Phyla","genomad_Class","genomad_Order","genomad_Family"))

# Check and split tax correctly through factor(levels)
for(i in 1:length(helmi$viral_name)){
  if(grepl("idae", helmi$genomad_Realm[i])){
    helmi$genomad_Family[i] <- helmi$genomad_Realm[i]
    helmi$genomad_Realm[i] <- NA
  }
  if(grepl("etes", helmi$genomad_Realm[i])){
    helmi$genomad_Class[i] <- helmi$genomad_Realm[i]
    helmi$genomad_Realm[i] <- NA
  }
  if(grepl("idae", helmi$genomad_Order[i])){
    helmi$genomad_Family[i] <- helmi$genomad_Order[i]
    helmi$genomad_Order[i] <- NA
  }
  if(grepl("ales", helmi$genomad_Kingdom[i])){
    helmi$genomad_Order[i] <- helmi$genomad_Kingdom[i]
    helmi$genomad_Kingdom[i] <- NA
  }
}

# Merge phagcn
phagcn <- phagcn %>% 
  rename(c("viral_name" = "contig_name","phagcn_Family" = "prediction")) %>%
  filter(viral_name %in% helmi$viral_name) %>%
  select(c("viral_name","phagcn_Family"))

helmi <- merge(helmi, phagcn, all=T)

# Keep only prediction in cuado class
helmi <- helmi %>% 
  mutate(
    phagcn_Family = case_when(
      (helmi$genomad_Class != "Caudoviricetes" | is.na(helmi$genomad_Class)) ~ NA,
      .default = helmi$phagcn_Family
    )
  )

# Phatyp
phatyp <- phatyp %>% 
  rename(c("viral_name" = "Contig","phatyp_pred" = "Pred"))

phatyp_sub <- subset(phatyp, viral_name %in% helmi$viral_name)
helmi <- merge(helmi, phatyp_sub, all=T)

helmi <- helmi %>% 
  mutate(
    lifestyle = case_when(
      ((helmi$phatyp_pred == "temperate") & (helmi$Score >= 0.7)) ~ "Temperate",
      ((helmi$phatyp_pred == "temperate") & (helmi$Score < 0.7) & ((helmi$provirus == "Yes") & (helmi$topology == "Provirus"))) ~ "Temperate",
      ((helmi$phatyp_pred == "temperate") & (helmi$Score < 0.7) &  
         (((helmi$provirus != "Yes") & (helmi$topology == "Provirus")) |
            ((helmi$provirus == "Yes") & (helmi$topology != "Provirus")) |
            ((helmi$provirus != "Yes") & (helmi$topology != "Provirus")))) ~ "Uncertain Temperate",
      ((helmi$phatyp_pred == "virulent") & (helmi$Score < 0.7)) ~ "Uncertain Virulent",
      ((helmi$phatyp_pred == "virulent") & (helmi$Score >= 0.7)) ~ "Virulent"
    )
  )

# Add lifestyle evidence 
helmi <- helmi %>% 
  mutate(
    lifestyle_evidence = case_when(
      (helmi$lifestyle == "Virulent" | helmi$lifestyle == "Uncertain Virulent") ~ "PhaTyp",
      ((helmi$phatyp_pred == "temperate") & (helmi$Score >= 0.7) & (helmi$provirus == "Yes") & (helmi$topology == "Provirus")) ~ "PhaTyp/CheckV/Genomad",
      (((helmi$phatyp_pred == "temperate") & (helmi$Score < 0.7)) | ((helmi$provirus == "Yes") & (helmi$topology == "Provirus"))) ~ "CheckV/Genomad",
      ((helmi$lifestyle == "Temperate") &  
         (((helmi$provirus != "Yes") & (helmi$topology == "Provirus")) |
            ((helmi$provirus == "Yes") & (helmi$topology != "Provirus")) | 
            ((helmi$provirus == "No") & (helmi$topology != "Provirus")))) ~ "PhaTyp",
      ((helmi$lifestyle == "Uncertain Temperate") &  
         (((helmi$provirus != "Yes") & (helmi$topology == "Provirus")) |
            ((helmi$provirus == "Yes") & (helmi$topology != "Provirus")) | 
            ((helmi$provirus == "No") & (helmi$topology != "Provirus")))) ~ "PhaTyp"
    )
  )

helmi_tools <- helmi %>% filter(DVF == T & VS2 == T)
helmi_cycle <- helmi %>% group_by(lifestyle) %>% tally()

# Remove non-phages / Ensure only phages are annotated - Caudoviricetes,Faserviricetes,Malgrandaviricetes
# Keep >5000 length
tax_count <- helmi %>% count(genomad_Class)
helmi <- helmi %>% 
  filter(genomad_Class == "Caudoviricetes" | genomad_Class == "Faserviricetes" | genomad_Class == "Malgrandaviricetes")
# Tectiliviricetes is both but at the Family level that taxa is eukarytoic
length_small <- helmi %>% filter(contig_length < 5000)
helmi <- helmi %>% filter(contig_length > 4999)
#write_xlsx(helmi, "helmi_catalog.xlsx", col_names=TRUE)

helmi_clust <- read_tsv("c95_species_vOTU_cluster.tsv", col_names = F)

# Samples
vOTUs <- helmi$viral_name
vOTU_clust <- subset(helmi_clust, X1 %in% vOTUs)

vOTU_grouped <- vOTU_clust %>%
  group_by(X1) %>%
  summarise(sample = toString(X2), .groups = 'drop') 

wanted_samples <- paste0("M|F|B4|B5|B7|B9")
vOTU_grouped <- vOTU_grouped %>% 
  mutate(sample = map(sample, ~unique(unlist(str_extract_all(.x, wanted_samples))))) 

helmi <- helmi %>% left_join(vOTU_grouped, by = c("viral_name" = "X1"))

# 
helmi <- helmi %>% 
  unnest(sample) %>%
  group_by(viral_name,contig_length,DVF,VS2,provirus,proviral_length,checkv_quality,miuvig_quality,warnings,topology,genomad_Virus,genomad_Realm,
           genomad_Kingdom,genomad_Phyla,genomad_Class,genomad_Order,genomad_Family,phagcn_Family,phatyp_pred,Score,lifestyle,lifestyle_evidence) %>% 
  summarise(sample = paste(sample, collapse = ",")) 

helmi <- helmi[,c(1,23,2:22)]
#write_xlsx(helmi, "helmi_catalog_fin.xlsx", col_names=TRUE)

helmi <- read_excel("helmi_catalog_fin.xlsx")

# CheckV quality 
checkv_count <- helmi %>% count(checkv_quality)

# Taxonomy
tax <- helmi %>% count(genomad_Class)

# lifestyle 
lifestlye <- helmi %>% count(lifestyle)

# Plots 

Pal6 <- c("#BB3E03", "#E9C46A", "#CA6702", "#0A9396",
          "#005F73", "#80BBA8")

# plot quality
checkv_count$checkv_quality <- factor(checkv_count$checkv_quality, levels=c('Complete','High-quality','Medium-quality','Low-quality','Not-determined'))
quality <- ggplot(checkv_count, aes(x="", y=n, fill=checkv_quality))+
  geom_bar(width = 1, stat = "identity") + theme_classic() + 
  labs(x = "", y = "Number Sequence") +
  ggtitle("Quality of Sequences in HVrC") +
  scale_fill_manual(values = Pal6)  +
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 10),legend.title = element_text(size = 10), axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))
quality
#ggsave("quality.pdf", width = 5, height = 4)

colors <- distinctColorPalette(3)
# plot taxonomy
taxonomy <- ggplot(tax, aes(x = 1, y = n, fill = genomad_Class)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme(panel.background = element_blank(), plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()) +
  ggtitle("HELMI Library Viral Class") + 
  scale_fill_manual(values = colors, guide = "none") +
  guides(fill = guide_legend(title = "Class"))
taxonomy
#ggsave("class_chart.pdf")

lifestyle_plot <- helmi %>% group_by(lifestyle)  %>% tally()

# plot lifestyle
cycle <- ggplot(lifestyle_plot, aes(x = 1, y = n, fill = lifestyle)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme(panel.background = element_blank(), plot.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle("HELMI Library Viral Lifestyle") + 
  scale_fill_manual(values = Pal6, guide = "none") +
  guides(fill = guide_legend(title = "Class"))
cycle
ggsave("lifestyle_chart.pdf")

tax2 <- tax %>% mutate(Percentage = (n/sum(n))*100) %>%
  select(c("genomad_Class","Percentage")) %>%
  rename("Class" = "genomad_Class")
  
#pdf("Class_table.pdf", height=4, width=4)
#grid.table(tax)
#dev.off()
  
# DVF check 
DVF_test <- genomad_missing %>% 
  mutate(
    ID = case_when(
      (DVF == T) & (VS2 == T) ~ "Both",
      (DVF == T) & (VS2 == F) ~ "DVF",
      (DVF == F) & (VS2 == T) ~ "VS2",
      (DVF == F) & (VS2 == F) ~ "None"
    )
  )

DVF_test_count <- DVF_test %>% group_by(ID) %>% count()
DVF_test_2 <- DVF_test %>% filter(ID == "Both" & contig_length > 4999)
