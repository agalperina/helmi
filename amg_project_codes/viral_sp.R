library(ape)
library(dplyr)
library(phyloseq)
library(tidyverse)
library(microbiome)
library(readxl)
library(ggpubr)
library(stringi)
library(vegan)

setwd("~/OneDrive/University of Helsinki/msc_thesis/viral_species")
helmi_vr <- read_xlsx("../viral_selection_2/helmi_catalog_fin.xlsx")

order_pal <- c("#3890BC", "#1C489A", "#000077","#808080")

Pal6 <- c("#BB3E03", "#E9C46A", "#CA6702", "#0A9396",
          "#005F73", "#80BBA8")

Pal12 <- c("#9B2226","#CC4400", "#D66915", "#E08E29", "#F0C761", "#FFFF99",
           "#C2FCFF", "#7CC6DE", "#3890BC", "#1C489A", "#000077", "#001219")

Pal20 <- c("#FFBA08","#FAA307","#F48C06","#E85D04","#DC2F02","#D00000","#9D0208","#6A040F","#370617","#03071E","#012A4A","#013A63", "#01497C", "#014F86", "#2A6F97", "#2C7DA0",
           "#468FAF", "#61A5C2", "#89C2D9", "#808080")
Pal_20_2 <- c("#FFBA08","#89C2D9","#FAA307", "#61A5C2","#F48C06","#468FAF","#E85D04","#2C7DA0","#DC2F02","#2A6F97","#D00000","#014F86","#9D0208", "#01497C","#6A040F","#013A63","#370617","#012A4A","#03071E","#808080")

# generate phyloseq obj
# subset catalog for tax table 
taxonomy <- helmi_vr %>%
  select(viral_name, genomad_Realm:genomad_Class, phagcn_Family) %>%
  column_to_rownames("viral_name") %>%
  rename("Realm" = "genomad_Realm", "Kingdom" = "genomad_Kingdom","Phyla" = "genomad_Phyla","Class" = "genomad_Class","Family" = "phagcn_Family") 

TAX = tax_table(as.matrix(taxonomy))

taxonomy_o <- helmi_vr %>%
  select(viral_name, genomad_Realm:genomad_Order, phagcn_Family) %>%
  column_to_rownames("viral_name") %>%
  rename("Realm" = "genomad_Realm", "Kingdom" = "genomad_Kingdom","Phyla" = "genomad_Phyla","Class" = "genomad_Class","Order" = "genomad_Order","Family" = "phagcn_Family") 
    
TAX_o = tax_table(as.matrix(taxonomy_o)) # also create order level physeq obj to laster look at crassvirales viruses

# create count table
otu_table <- read_csv("All_Sec3Counts_table.csv")
otuMat <- otu_table %>% tibble::column_to_rownames("Name")
OTU = otu_table(as.matrix(otuMat), taxa_are_rows = TRUE)

# sample data
sampleMat <- read_csv("sample_desc_table_01.05.22.csv") 
sampleMat <- sampleMat %>%tibble::column_to_rownames("Sample_ID")
SAMPLE = sample_data(as.data.frame(sampleMat))

physeq = phyloseq(OTU, TAX, SAMPLE)
physeq_o = phyloseq(OTU, TAX_o, SAMPLE)
#saveRDS(physeq, "PhySeq_helmi_viral_27.11.23.rds")

##
physeq <- readRDS("PhySeq_helmi_viral_27.11.23.rds")

# Alpha Diversity 

#tab <- microbiome::alpha(physeq, index = "all") 
#saveRDS(tab, "viral_alpha_diversity.rds")
tab <- readRDS("viral_alpha_diversity.rds")
tab_sep <- tab %>% 
  rownames_to_column(var="Sample_ID") %>%
  separate(Sample_ID, into=c("Family_ID", "Sample_type"), sep="_", remove=FALSE)
tab_sep$Sample_type <- gsub('M', 'P',
                            gsub('F', 'P', tab_sep$Sample_type))

comp <- list(c("B4","B5"),c("B5","B7"), c("B7","B9"),c("B9","P"))
alpha_plot <- tab_sep %>% ggplot(aes(x=Sample_type, y=diversity_shannon, color=Sample_type))+
  geom_boxplot()+
  theme_classic()+
  scale_color_manual(values=Pal6) +
  theme(legend.position = "none") + 
  stat_compare_means(comparisons = comp, method = "wilcox.test")
alpha_plot
ggsave("viral_alpha.pdf")

# Richness
rich_plot <- tab_sep %>% ggplot(aes(x=Sample_type, y=chao1, color=Sample_type))+
  geom_boxplot()+
  theme_light()+
  scale_color_manual(values=Pal6) +
  theme(legend.position = "none") +
  theme_classic()+
  stat_compare_means(comparisons = comp, method = "wilcox.test")
rich_plot
ggsave("viral_rich.pdf")

# Avrg composition at family level
cluster_B4.5 <- read_csv("../bact_comm/HClust_F_CLR_B4.5_18.08.22.csv")
cluster_B4.5$ClusterID <- sub("^", "EC_", cluster_B4.5$ClusterID)
cluster_B7.9 <- read_csv("../bact_comm/HClust_F_CLR_B7.9_18.08.22.csv")
cluster_B7.9$ClusterID <- sub("^", "LC_", cluster_B7.9$ClusterID)
cluster_P <- read_csv("../bact_comm/HClust_F_CLR_P_13.11.23.csv") %>% rename('ClusterID' = "ClusterID_allWGS")
cluster_P$ClusterID <- sub("^", "P_", cluster_P$ClusterID)

cluster_info <- rbind(cluster_B4.5,cluster_B7.9,cluster_P)
# relative abundance composition plots
physeq_comp  <- microbiome::transform(physeq, "compositional") 
#psG  <- aggregate_rare(physeq_comp, level="Phyla", detection=0.1, prevalence=0.1) 
physeq_fam <- aggregate_taxa(physeq_comp, level="Family") # should I do aggregate_rar here?
plot_data <- psmelt(physeq_fam) %>% select(Sample, Sample_type, OTU, Abundance)
plot_data <- left_join(plot_data,cluster_info, by = c("Sample" = "Sample_ID"))

plot_data$Sample_type <- gsub('M', 'P', gsub('F', 'P', plot_data$Sample_type))

plot_data_av <- plot_data %>% group_by(OTU, Sample_type,ClusterID) %>% 
  mutate(av_Abundance=mean(Abundance)) %>%
  select(OTU, Sample_type, av_Abundance,ClusterID) %>% unique()

plot_data_av$OTU <- str_replace(plot_data_av$OTU, "Unknown", "zzzUnclassified")
other <- c("Zobellviridae", "Zierdtviridae", "Schitoviridae", "Straboviridae", 
           "Routreeviridae", "Kyanoviridae", "Demerecviridae", "Chaseviridae", "Autographiviridae")
plot_data_av$OTU <- stri_replace_all_regex(plot_data_av$OTU ,
                                  pattern=other,
                                  replacement="zzOther",
                                  vectorize=FALSE)
pal13 <- c("#9B2226","#AE2012", "#BB3E03", "#CA6702", "#EE9B00", "#E9C46A",
           "#80BBA8", "#0A9396", "#005F73", "#033549","#012A4A", "#808080","#3d3d3d")
comp_plot_cluster <- plot_data_av %>%  ggplot(aes(fill=OTU, y=av_Abundance, x=ClusterID)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = pal13) +
  theme_bw() +
  facet_grid(.~Sample_type, scales = "free_x")
comp_plot_cluster
ggsave("viral_comp_plot_cluster_13.2.24.pdf")

# remove cluster info 
plot_data <- psmelt(physeq_fam) %>% select(Sample, Sample_type, OTU, Abundance)
plot_data$Sample_type <- gsub('M', 'P', gsub('F', 'P', plot_data$Sample_type))

plot_data_av <- plot_data %>% group_by(OTU, Sample_type) %>% 
  mutate(av_Abundance=mean(Abundance)) %>%
  select(OTU, Sample_type, av_Abundance) %>% unique()

plot_data_av$OTU <- str_replace(plot_data_av$OTU, "Unknown", "zzzUnclassified")
plot_data_av$OTU <- stri_replace_all_regex(plot_data_av$OTU ,
                                           pattern=other,
                                           replacement="zzOther",
                                           vectorize=FALSE)
comp_plot <- plot_data_av %>%  ggplot(aes(fill=OTU, y=av_Abundance, x=Sample_type)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = pal13) +
  theme_classic() 
comp_plot
dev.copy2pdf(file="viral_rel_comp_13.2.24.pdf")
ggsave("viral_rel_comp_13.2.24.pdf")

physeq_comp_o  <- microbiome::transform(physeq_o, "compositional") 
#psG  <- aggregate_rare(physeq_comp, level="Phyla", detection=0.1, prevalence=0.1) 
physeq_fam_o <- aggregate_taxa(physeq_comp_o, level="Order") 
plot_data_o <- psmelt(physeq_fam_o) %>% select(Sample, Sample_type, OTU, Abundance)

plot_data_av_o <- plot_data_o %>% group_by(OTU, Sample_type) %>% 
  mutate(av_Abundance=mean(Abundance)) %>%
  select(OTU, Sample_type, av_Abundance) %>% unique() %>%
  filter(OTU != "Unknown")

plot_data_av_o$OTU <- str_replace(plot_data_av_o$OTU, "Unknown", "Unclassified")

## Crass phages
physeq_crass <- aggregate_taxa(physeq_o, level="Order") 
physeq_data_crass <- psmelt(physeq_crass) %>% select(Sample, Sample_type, OTU, Abundance) %>%  filter(OTU == "Crassvirales")

physeq_data_crass$Sample_type <- gsub('M', 'Parental',
                             gsub('F', 'Parental', physeq_data_crass$Sample_type))

b4 <- helmi_vr %>% filter(grepl("B4",viral_name)) %>% filter(genomad_Order == "Crassvirales")
b9 <- helmi_vr %>% filter(grepl("B9",viral_name)) %>% filter(genomad_Order == "Crassvirales")
m <- helmi_vr %>% filter(grepl("M",viral_name)) %>% filter(genomad_Order == "Crassvirales")
f <- helmi_vr %>% filter(grepl("F",viral_name)) %>% filter(genomad_Order == "Crassvirales")


comp <- list(c("B4","B5"),c("B5","B7"),c("B7","B9"),c("B9","Parental"))
crass_plot <- physeq_data_crass %>%  ggplot(aes(color=Sample_type, y=log(Abundance+0.0001), x=Sample_type)) + 
  geom_boxplot() + scale_color_manual(values = Pal6) +
  theme_classic() +
  stat_compare_means(comparisons = comp, method = "wilcox.test")
crass_plot
ggsave("crass_ab.pdf")

physeq_data_crass_b4 <- physeq_data_crass %>% filter(Sample_type == "B4")
# PcoA clr + euclidean (Aitchinson)

# Data prep
taxa <- core_members(physeq_comp, detection = 0.1/100, prevalence = 1/100)
taxa_unique <- unique(taxa)
physeq_denoised <- prune_taxa(taxa_unique, physeq)
ps_clr_abund <- microbiome::transform(physeq_denoised, "clr")

# Aitchinson computation
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_clr_abund))
ps_rel_otu <- t(ps_rel_otu)
ait_dist <- vegan::vegdist(ps_rel_otu, method = "euclidean")

ord <- ordinate(physeq_denoised, "MDS", distance = "ait_dist")
#ord_bray <- ordinate(physeq_denoised, "MDS", "bray")
Pal5 <- c("#9B2226", "#BB3E03", "#EE9B00", "#E9D8A6",
          "#005F73")
physeq_denoised@sam_data$Sample_type <- gsub('M', 'P', gsub('F', 'P', physeq_denoised@sam_data$Sample_type))
# plot 
plot_ordination(physeq_denoised, ord, color = "Sample_type") + 
  geom_point(size = 0.75) +  
  scale_color_manual(values = Pal5) + 
  theme_classic()
ggsave("viral_PCoA_all_timepoints.pdf")

perm_metadata <- data.frame(sample_data(physeq_denoised))
disp <- anova(betadisper(ait_dist, perm_metadata$Sample_type))
disp

## Lifestyle
helmi_sub_l <- helmi_vr %>% select(viral_name,lifestyle)
helmi_sub_l$lifestyle <- helmi_sub_l$lifestyle %>% replace_na('Undetermined')
F_count <- read_csv2("F_sel3count_table.csv")
F_lifestyle <- left_join(F_count, helmi_sub_l, by = c("Name" = "viral_name")) 
F_lifestyle <- F_lifestyle[,c(1,125,2:124)]

# sum the samples per lifestyle
F_lifestyle_l <- F_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

F_lifestly_avg <- F_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

F_lifestly_avg <- aggregate( .~lifestyle, data = F_lifestly_avg, FUN = mean)
#
M_count <- read_csv2("M_sel3count_table.csv")
M_lifestyle <- left_join(M_count, helmi_sub_l, by = c("Name" = "viral_name")) 
M_lifestyle <- M_lifestyle[,c(1,306,2:305)]

# sum the samples per lifestyle
M_lifestyle_l <- M_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

M_lifestly_avg <- M_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

M_lifestly_avg <- aggregate( .~lifestyle, data = M_lifestly_avg, FUN = mean)
#
B4_count <- read_csv2("B4_sel3count_table.csv")
B4_lifestyle <- left_join(B4_count, helmi_sub_l, by = c("Name" = "viral_name")) 
B4_lifestyle <- B4_lifestyle[,c(1,477,2:476)]

# sum the samples per lifestyle
B4_lifestyle_l <- B4_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B4_lifestly_avg <- B4_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

B4_lifestly_avg <- aggregate( .~lifestyle, data = B4_lifestly_avg, FUN = mean)
B4_lifestly_avg <- B4_lifestly_avg %>% mutate(prop = sum/(sum(B4_lifestly_avg$sum)))
#
B5_count <- read_csv2("B5_sel3count_table.csv")
B5_lifestyle <- left_join(B5_count, helmi_sub_l, by = c("Name" = "viral_name")) 
B5_lifestyle <- B5_lifestyle[,c(1,477,2:476)]

# sum the samples per lifestyle
B5_lifestyle_l <- B5_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B5_lifestly_avg <- B5_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

B5_lifestly_avg <- aggregate( .~lifestyle, data = B5_lifestly_avg, FUN = mean)
B5_lifestly_avg <- B5_lifestly_avg %>% mutate(prop = sum/(sum(B5_lifestly_avg$sum)))
#
B7_count <- read_csv2("B7_sel3count_table.csv")
B7_lifestyle <- left_join(B7_count, helmi_sub_l, by = c("Name" = "viral_name")) 
B7_lifestyle <- B7_lifestyle[,c(1,476,2:475)]

# sum the samples per lifestyle
B7_lifestyle_l <- B7_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B7_lifestly_avg <- B7_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

B7_lifestly_avg <- aggregate( .~lifestyle, data = B7_lifestly_avg, FUN = mean)
B7_lifestly_avg <- B7_lifestly_avg %>% mutate(prop = sum/(sum(B7_lifestly_avg$sum)))

#
B9_count <- read_csv2("B9_sel3count_table.csv")
B9_lifestyle <- left_join(B9_count, helmi_sub_l, by = c("Name" = "viral_name")) 
B9_lifestyle <- B9_lifestyle[,c(1,477,2:476)]

# sum the samples per lifestyle
B9_lifestyle_l <- B9_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B9_lifestly_avg <- B9_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

B9_lifestly_avg <- aggregate( .~lifestyle, data = B9_lifestly_avg, FUN = mean)
B9_lifestly_avg <- B9_lifestly_avg %>% mutate(prop = sum/(sum(B9_lifestly_avg$sum)))

# plotting data (combine)
P_lifestly_avg <- rbind(F_lifestly_avg,M_lifestly_avg)
P_lifestly_avg <- P_lifestly_avg %>% mutate(prop = sum/(sum(P_lifestly_avg$sum)))

lifestyle_plot <- rbind(B4_lifestly_avg,B5_lifestly_avg,B7_lifestly_avg,B9_lifestly_avg,P_lifestly_avg)
lifestyle_plot$sample <- rep(c("B4","B5","B7","B9","P"),times=c(5,5,5,5,10))

lifestyle_plot$lifestyle <- str_replace(lifestyle_plot$lifestyle, "Undetermined", "zzUndetermined")

Pal5 <- c("#BB3E03", "#CA6702", "#0A9396",
          "#005F73", "#808080")
lifestyle <- lifestyle_plot %>%  ggplot(aes(fill=lifestyle, y=prop, x=sample)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = Pal5) +
  theme_classic() 
lifestyle
ggsave("lifestyle_2_13_24.pdf")

## Temperate rel abd
test <- sweep(B9_count[-1], 2, colSums(B9_count[,-1]), '/') * 100
virals_name <- B9_count[1]
B9_count <- cbind(virals_name, test)
B9_lifestyle <- left_join(B9_count, helmi_sub_l, by = c("Name" = "viral_name")) 
B9_lifestyle <- B9_lifestyle[,c(1,477,2:476)]

# sum the samples per lifestyle
B9_lifestyle_l <- B9_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 
B9_temp_avg <- B9_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample) %>%
  filter(grepl("Virulent",lifestyle))

test <- sweep(B7_count[-1], 2, colSums(B7_count[,-1]), '/') * 100
virals_name <- B7_count[1]
B7_count <- cbind(virals_name, test)
B7_lifestyle <- left_join(B7_count, helmi_sub_l, by = c("Name" = "viral_name")) 
B7_lifestyle <- B7_lifestyle[,c(1,476,2:475)]

# sum the samples per lifestyle
B7_lifestyle_l <- B7_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B7_temp_avg <- B7_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample) %>%
  filter(grepl("Virulent",lifestyle))

test <- sweep(B5_count[-1], 2, colSums(B5_count[,-1]), '/') * 100
virals_name <- B5_count[1]
B5_count <- cbind(virals_name, test)
B5_lifestyle <- left_join(B5_count, helmi_sub_l, by = c("Name" = "viral_name")) 
B5_lifestyle <- B5_lifestyle[,c(1,477,2:476)]

# sum the samples per lifestyle
B5_lifestyle_l <- B5_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B5_temp_avg <- B5_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample) %>%
  filter(grepl("Virulent",lifestyle))

test <- sweep(B4_count[-1], 2, colSums(B4_count[,-1]), '/') * 100
virals_name <- B4_count[1]
B4_count <- cbind(virals_name, test)
B4_lifestyle <- left_join(B4_count, helmi_sub_l, by = c("Name" = "viral_name")) 
B4_lifestyle <- B4_lifestyle[,c(1,477,2:476)]

# sum the samples per lifestyle
B4_lifestyle_l <- B4_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B4_temp_avg <- B4_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample) %>%
  filter(grepl("Virulent",lifestyle))

test <- sweep(M_count[-1], 2, colSums(M_count[,-1]), `/`) * 100
virals_name <- M_count[1]
M_count <- cbind(virals_name, test)
M_lifestyle <- left_join(M_count, helmi_sub_l, by = c("Name" = "viral_name")) 
M_lifestyle <- M_lifestyle[,c(1,306,2:305)]

# sum the samples per lifestyle
M_lifestyle_l <- M_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

M_temp_avg <- M_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample) %>%
  filter(grepl("Virulent",lifestyle))

test <- sweep(F_count[-1], 2, colSums(F_count[,-1]), `/`) * 100
virals_name <- F_count[1]
F_count <- cbind(virals_name, test)
F_lifestyle <- left_join(F_count, helmi_sub_l, by = c("Name" = "viral_name")) 
F_lifestyle <- F_lifestyle[,c(1,125,2:124)]

# sum the samples per lifestyle
F_lifestyle_l <- F_lifestyle %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

F_temp_avg <- F_lifestyle_l %>%
  group_by(lifestyle, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample) %>%
  filter(grepl("Virulent",lifestyle))

# merge parents
temp_lifestyle_plot <- rbind(B4_temp_avg,B5_temp_avg,B7_temp_avg,B9_temp_avg,M_temp_avg,F_temp_avg)
temp_lifestyle_plot$sample <- rep(c("B4","B5","B7","B9","M","F"),times=c(950,950,948,950,608,246))

comp_1 <- list(c("B4","B5"),c("B5","B7"),c("B7","B9"),c("B9","M"),c("B9","F"))
comp_2 <- list(c("B9","F"),c("B9","M"),c("B7","F"),c("B7","M"),c("B5","F"),c("B5","M"),c("B4","F"),c("B4","M"))
# plot data
temp_avg_plot <- temp_lifestyle_plot %>%  ggplot(aes(color=sample, y=sum, x=sample)) + 
  geom_boxplot() + scale_color_manual(values = Pal6) +
  theme_classic() +
  stat_compare_means(comparisons = comp_1, method = "wilcox.test")
temp_avg_plot



