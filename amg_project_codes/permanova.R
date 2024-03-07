library(ape)
library(dplyr)
library(phyloseq)
library(tidyverse)
library(microbiome)
library(readxl)
library(ggpubr)
library(gridExtra)
library(vegan)

# permanova
setwd("~/OneDrive/University of Helsinki/msc_thesis/viral_species")
helmi_vr <- read_xlsx("../viral_selection_2/helmi_catalog_fin.xlsx")

Pal3 <- c("#BB3E03", "#005F73", "#E9C46A")
Pal20 <- c("#FFBA08","#FAA307","#F48C06","#E85D04","#DC2F02","#D00000","#9D0208","#6A040F","#370617","#03071E","#012A4A","#013A63", "#01497C", "#014F86", "#2A6F97", "#2C7DA0",
           "#468FAF", "#61A5C2", "#89C2D9","#808080")

# read in FCT cluster info and combine
sampleMat <- read_csv("sample_desc_table_01.05.22.csv") 
cluster_B4.5 <- read_csv("../bact_comm/HClust_F_CLR_B4.5_18.08.22.csv")
cluster_B4.5$ClusterID <- sub("^", "EC_", cluster_B4.5$ClusterID)
cluster_B7.9 <- read_csv("../bact_comm/HClust_F_CLR_B7.9_18.08.22.csv")
cluster_B7.9$ClusterID <- sub("^", "LC_", cluster_B7.9$ClusterID)
cluster_P <- read_csv("../bact_comm/HClust_F_CLR_P_13.11.23.csv") %>% rename('ClusterID' = "ClusterID_allWGS")
cluster_P$ClusterID <- sub("^", "P_", cluster_P$ClusterID)

clusters <- rbind(cluster_B4.5,cluster_B7.9,cluster_P)

sampleMat <- left_join(clusters, sampleMat, by = "Sample_ID")

taxonomy <- helmi_vr %>%
  select(viral_name, genomad_Realm:genomad_Class, phagcn_Family) %>%
  column_to_rownames("viral_name") %>%
  rename("Realm" = "genomad_Realm", "Kingdom" = "genomad_Kingdom","Phyla" = "genomad_Phyla","Class" = "genomad_Class","Family" = "phagcn_Family") 

TAX = tax_table(as.matrix(taxonomy))

# create count table
otu_table <- read_csv("All_Sec3Counts_table.csv")
otuMat <- otu_table %>% tibble::column_to_rownames("Name")
OTU = otu_table(as.matrix(otuMat), taxa_are_rows = TRUE)

# sample data
sampleMat <- sampleMat %>%tibble::column_to_rownames("Sample_ID")
SAMPLE = sample_data(as.data.frame(sampleMat))

physeq = phyloseq(OTU, TAX, SAMPLE)

# Split by timepoint
physeq_B4 <- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "B4")
physeq_B5 <- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "B5")
physeq_B7 <- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "B7")
physeq_B9 <- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "B9")

physeq_B4_comp  <- microbiome::transform(physeq_B4, "compositional") 
physeq_B5_comp  <- microbiome::transform(physeq_B5, "compositional") 
physeq_B7_comp  <- microbiome::transform(physeq_B7, "compositional") 
physeq_B9_comp  <- microbiome::transform(physeq_B9, "compositional") 

# B4 
taxa_B4 <- core_members(physeq_B4_comp, detection = 0.1/100, prevalence = 1/100)
taxa_B4.unique <- unique(taxa_B4)
physeq_B4_denoised <- prune_taxa(taxa_B4.unique, physeq_B4)
ps_clr_abund_B4 <- microbiome::transform(physeq_B4_denoised, "clr")

ps_rel_otu_B4 <- data.frame(phyloseq::otu_table(ps_clr_abund_B4))
ps_rel_otu_B4 <- t(ps_rel_otu_B4)
ait_dist_B4 <- vegan::vegdist(ps_rel_otu_B4, method = "euclidean")

ord_B4 <- ordinate(physeq_B4_denoised, "PCoA", distance = "ait_dist_B4")

PCoA_B4 <- plot_ordination(physeq_B4_denoised, ord_B4, color = "ClusterID") + 
  theme_classic() +
  theme(strip.background = element_blank()) + scale_color_manual(values = Pal3)
PCoA_B4
#ggsave("viral_PCoA_B4.pdf")

B4_perm_metadata <- data.frame(sample_data(physeq_B4_denoised))
B4_adonis <- adonis2(ait_dist_B4 ~ ClusterID, data = B4_perm_metadata)
B4_adonis

#pdf("B4_adonis.pdf", height=8, width=8)
#grid.table(B4_adonis)
#dev.off()

B4_disp <- anova(betadisper(ait_dist_B4, B4_perm_metadata$ClusterID))
B4_disp

#pdf("B4_beta_disp.pdf", height=10, width=10)
#grid.table(B4_disp)
#dev.off()

# B5
taxa_B5 <- core_members(physeq_B5_comp, detection = 0.1/100, prevalence = 1/100)
taxa_B5.unique <- unique(taxa_B5)
physeq_B5_denoised <- prune_taxa(taxa_B5.unique, physeq_B5)
ps_clr_abund_B5 <- microbiome::transform(physeq_B5_denoised, "clr")

ps_rel_otu_B5 <- data.frame(phyloseq::otu_table(ps_clr_abund_B5))
ps_rel_otu_B5 <- t(ps_rel_otu_B5)
ait_dist_B5 <- vegan::vegdist(ps_rel_otu_B5, method = "euclidean")

ord_B5 <- ordinate(physeq_B5_denoised, "PCoA", distance = "ait_dist_B5")

PCoA_B5 <- plot_ordination(physeq_B5_denoised, ord_B5, color = "ClusterID") + 
  theme_classic() +
  theme(strip.background = element_blank()) + scale_color_manual(values = Pal3)
PCoA_B5
ggsave("viral_PCoA_B5.pdf")

B5_perm_metadata <- data.frame(sample_data(physeq_B5_denoised))
B5_adonis <- adonis2(ait_dist_B5 ~ ClusterID, data = B5_perm_metadata)
B5_adonis

pdf("B5_adonis.pdf", height=8, width=8)
grid.table(B5_adonis)
dev.off()

B5_disp <- anova(betadisper(ait_dist_B5, B5_perm_metadata$ClusterID))
B5_disp

pdf("B5_beta_disp.pdf", height=10, width=10)
grid.table(B5_disp)
dev.off()

# B7
taxa_B7 <- core_members(physeq_B7_comp, detection = 0.1/100, prevalence = 1/100)
taxa_B7.unique <- unique(taxa_B7)
physeq_B7_denoised <- prune_taxa(taxa_B7.unique, physeq_B7)
ps_clr_abund_B7 <- microbiome::transform(physeq_B7_denoised, "clr")

ps_rel_otu_B7 <- data.frame(phyloseq::otu_table(ps_clr_abund_B7))
ps_rel_otu_B7 <- t(ps_rel_otu_B7)
ait_dist_B7 <- vegan::vegdist(ps_rel_otu_B7, method = "euclidean")

ord_B7 <- ordinate(physeq_B7_denoised, "PCoA", distance = "ait_dist_B7")

PCoA_B7 <- plot_ordination(physeq_B7_denoised, ord_B7, color = "ClusterID") + 
  theme_classic() +
  theme(strip.background = element_blank()) + scale_color_manual(values = Pal3)
PCoA_B7
ggsave("viral_PCoA_B7.pdf")

B7_perm_metadata <- data.frame(sample_data(physeq_B7_denoised))
B7_adonis <- adonis2(ait_dist_B7 ~ ClusterID, data = B7_perm_metadata)
B7_adonis

pdf("B7_adonis.pdf", height=8, width=8)
grid.table(B7_adonis)
dev.off()

B7_disp <- anova(betadisper(ait_dist_B7, B7_perm_metadata$ClusterID))
B7_disp

pdf("B7_beta_disp.pdf", height=10, width=10)
grid.table(B7_disp)
dev.off()

# B9
taxa_B9 <- core_members(physeq_B9_comp, detection = 0.1/100, prevalence = 1/100)
taxa_B9.unique <- unique(taxa_B9)
physeq_B9_denoised <- prune_taxa(taxa_B9.unique, physeq_B9)
ps_clr_abund_B9 <- microbiome::transform(physeq_B9_denoised, "clr")

ps_rel_otu_B9 <- data.frame(phyloseq::otu_table(ps_clr_abund_B9))
ps_rel_otu_B9 <- t(ps_rel_otu_B9)
ait_dist_B9 <- vegan::vegdist(ps_rel_otu_B9, method = "euclidean")

ord_B9 <- ordinate(physeq_B9_denoised, "PCoA", distance = "ait_dist_B9")

PCoA_B9 <- plot_ordination(physeq_B9_denoised, ord_B9, color = "ClusterID") + 
  theme_classic() +
  theme(strip.background = element_blank()) + scale_color_manual(values = Pal3)
PCoA_B9
ggsave("viral_PCoA_B9.pdf")

B9_perm_metadata <- data.frame(sample_data(physeq_B9_denoised))
B9_adonis <- adonis2(ait_dist_B9 ~ ClusterID, data = B9_perm_metadata)
B9_adonis

pdf("B9_adonis.pdf", height=8, width=8)
grid.table(B9_adonis)
dev.off()

B9_disp <- anova(betadisper(ait_dist_B9, B9_perm_metadata$ClusterID))
B9_disp

pdf("B9_beta_disp.pdf", height=10, width=10)
grid.table(B9_disp)
dev.off()

# At Family Level

physeq_B4_fam <- aggregate_taxa(physeq_B4_comp, level="Family") 
physeq_B5_fam <- aggregate_taxa(physeq_B5_comp, level="Family") 
physeq_B7_fam <- aggregate_taxa(physeq_B7_comp, level="Family") 
physeq_B9_fam <- aggregate_taxa(physeq_B9_comp, level="Family") 
   
# B4
taxa_B4_fam <- core_members(physeq_B4_fam, detection = 0.1/100, prevalence = 1/100)
taxa_B4.unique_fam <- unique(taxa_B4_fam)
physeq_B4_denoised_fam <- prune_taxa(taxa_B4.unique_fam, physeq_B4_fam)
ps_clr_abund_B4_fam <- microbiome::transform(physeq_B4_denoised_fam, "clr")

ps_rel_otu_B4_fam <- data.frame(phyloseq::otu_table(ps_clr_abund_B4_fam))
ps_rel_otu_B4_fam <- t(ps_rel_otu_B4_fam)
ait_dist_B4_fam <- vegan::vegdist(ps_rel_otu_B4_fam, method = "euclidean")

ord_B4_fam <- ordinate(physeq_B4_denoised_fam, "PCoA", distance = "ait_dist_B4_fam")

PCoA_B4_fam <- plot_ordination(physeq_B4_denoised_fam, ord_B4_fam, color = "ClusterID") + 
  theme_classic() +
  theme(strip.background = element_blank()) + scale_color_manual(values = Pal3)
PCoA_B4_fam
ggsave("viral_PCoA_B4_fam.pdf")

B4_perm_metadata_fam <- data.frame(sample_data(physeq_B4_denoised_fam))
B4_adonis_fam <- adonis2(ait_dist_B4_fam ~ ClusterID, data = B4_perm_metadata_fam)
B4_adonis_fam

pdf("B4_adonis_fam.pdf", height=8, width=8)
grid.table(B4_adonis_fam)
dev.off()

B4_disp_fam <- anova(betadisper(ait_dist_B4_fam, B4_perm_metadata_fam$ClusterID))
B4_disp_fam

pdf("B4_disp_fam.pdf", height=10, width=10)
grid.table(B4_disp_fam)
dev.off()

# B5
taxa_B5_fam <- core_members(physeq_B5_fam, detection = 0.1/100, prevalence = 1/100)
taxa_B5.unique_fam <- unique(taxa_B5_fam)
physeq_B5_denoised_fam <- prune_taxa(taxa_B5.unique_fam, physeq_B5_fam)
ps_clr_abund_B5_fam <- microbiome::transform(physeq_B5_denoised_fam, "clr")

ps_rel_otu_B5_fam <- data.frame(phyloseq::otu_table(ps_clr_abund_B5_fam))
ps_rel_otu_B5_fam <- t(ps_rel_otu_B5_fam)
ait_dist_B5_fam <- vegan::vegdist(ps_rel_otu_B5_fam, method = "euclidean")

ord_B5_fam <- ordinate(physeq_B5_denoised_fam, "PCoA", distance = "ait_dist_B5_fam")

PCoA_B5_fam <- plot_ordination(physeq_B5_denoised_fam, ord_B5_fam, color = "ClusterID") + 
  theme_classic() +
  theme(strip.background = element_blank()) + scale_color_manual(values = Pal3)
PCoA_B5_fam
ggsave("viral_PCoA_B5_fam.pdf")

B5_perm_metadata_fam <- data.frame(sample_data(physeq_B5_denoised_fam))
B5_adonis_fam <- adonis2(ait_dist_B5_fam ~ ClusterID, data = B5_perm_metadata_fam)
B5_adonis_fam

pdf("B5_adonis_fam.pdf", height=8, width=8)
grid.table(B5_adonis_fam)
dev.off()

B5_disp_fam <- anova(betadisper(ait_dist_B5_fam, B5_perm_metadata_fam$ClusterID))
B5_disp_fam

pdf("B5_disp_fam.pdf", height=10, width=10)
grid.table(B5_disp_fam)
dev.off()

# B7
taxa_B7_fam <- core_members(physeq_B7_fam, detection = 0.1/100, prevalence = 1/100)
taxa_B7.unique_fam <- unique(taxa_B7_fam)
physeq_B7_denoised_fam <- prune_taxa(taxa_B7.unique_fam, physeq_B7_fam)
ps_clr_abund_B7_fam <- microbiome::transform(physeq_B7_denoised_fam, "clr")

ps_rel_otu_B7_fam <- data.frame(phyloseq::otu_table(ps_clr_abund_B7_fam))
ps_rel_otu_B7_fam <- t(ps_rel_otu_B7_fam)
ait_dist_B7_fam <- vegan::vegdist(ps_rel_otu_B7_fam, method = "euclidean")

ord_B7_fam <- ordinate(physeq_B7_denoised_fam, "PCoA", distance = "ait_dist_B7_fam")

PCoA_B7_fam <- plot_ordination(physeq_B7_denoised_fam, ord_B7_fam, color = "ClusterID") + 
  theme_classic() +
  theme(strip.background = element_blank()) + scale_color_manual(values = Pal3)
PCoA_B7_fam
ggsave("PCoA_B7_fam.pdf")

B7_perm_metadata_fam <- data.frame(sample_data(physeq_B7_denoised_fam))
B7_adonis_fam <- adonis2(ait_dist_B7_fam ~ ClusterID, data = B7_perm_metadata_fam)
B7_adonis_fam

pdf("B7_adonis_fam.pdf", height=8, width=8)
grid.table(B7_adonis_fam)
dev.off()

B7_disp_fam <- anova(betadisper(ait_dist_B7_fam, B7_perm_metadata_fam$ClusterID))
B7_disp_fam

pdf("B7_disp_fam.pdf", height=10, width=10)
grid.table(B7_disp_fam)
dev.off()

# B9
taxa_B9_fam <- core_members(physeq_B9_fam, detection = 0.1/100, prevalence = 1/100)
taxa_B9.unique_fam <- unique(taxa_B9_fam)
physeq_B9_denoised_fam <- prune_taxa(taxa_B9.unique_fam, physeq_B9_fam)
ps_clr_abund_B9_fam <- microbiome::transform(physeq_B9_denoised_fam, "clr")

ps_rel_otu_B9_fam <- data.frame(phyloseq::otu_table(ps_clr_abund_B9_fam))
ps_rel_otu_B9_fam <- t(ps_rel_otu_B9_fam)
ait_dist_B9_fam <- vegan::vegdist(ps_rel_otu_B9_fam, method = "euclidean")

ord_B9_fam <- ordinate(physeq_B9_denoised_fam, "PCoA", distance = "ait_dist_B9_fam")

PCoA_B9_fam <- plot_ordination(physeq_B9_denoised_fam, ord_B9_fam, color = "ClusterID") + 
  theme_classic() +
  theme(strip.background = element_blank()) + scale_color_manual(values = Pal3)
PCoA_B9_fam
ggsave("PCoA_B9_fam.pdf")

B9_perm_metadata_fam <- data.frame(sample_data(physeq_B9_denoised_fam))
B9_adonis_fam <- adonis2(ait_dist_B9_fam ~ ClusterID, data = B9_perm_metadata_fam)
B9_adonis_fam

pdf("B9_adonis_fam.pdf", height=8, width=8)
grid.table(B9_adonis_fam)
dev.off()

B9_disp_fam <- anova(betadisper(ait_dist_B9_fam, B9_perm_metadata_fam$ClusterID))
B9_disp_fam

pdf("B9_disp_fam.pdf", height=10, width=10)
grid.table(B9_disp_fam)
dev.off()

# Plot Family level conposition
# composition plots
physeq_comp  <- microbiome::transform(physeq, "compositional") %>%
  subset_samples(Sample_type %in% c("B4","B5","B7", "B9"))
psG  <- aggregate_taxa(physeq_comp, level="Family") 

plot_data <- psmelt(psG) %>% select(Sample, Sample_type, OTU, Abundance, ClusterID)

plot_data_av <- plot_data %>% group_by(OTU, Sample_type, ClusterID) %>% 
  mutate(av_Abundance=mean(Abundance)) %>%
  select(OTU, Sample_type, ClusterID, av_Abundance) %>% unique()

plot_data_av$OTU <- str_replace(plot_data_av$OTU, "Unknown", "zzzUnclassified")
plot_data_av$ClusterID <- gsub("post6month","LC",plot_data_av$ClusterID)
plot_data_av$ClusterID <- gsub("pre6month","EC",plot_data_av$ClusterID)

library(ggplot2)
library(plotly)

p <- plot_data_av %>%  ggplot(aes(fill=OTU, y=av_Abundance, x=ClusterID)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(.~Sample_type, scales = "free_x") +
  scale_fill_manual(values = Pal20) +
  ggtitle("Caudoviricetes Famile Level Timepoint Relative Abundance") + 
  theme_bw()
#ggsave("family_timepoint_abd.pdf")
ggplotly(p)

# checking to see how mnay samples are in each timepoint and determine missing samples
physeq_B4<- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "B4") # 475
physeq_B5 <- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "B5") # 475
physeq_B7<- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "B7") # 474
physeq_B9 <- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "B9") # 475
physeq_M <- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "M") # 304
physeq_F <- physeq %>%
  subset_samples(physeq@sam_data$Sample_type == "F") # 123

B4_samples <- physeq_B4@sam_data$Family_ID
B7_samples <- physeq_B7@sam_data$Family_ID #perhe1138

!(B4_samples %in% B7_samples)


