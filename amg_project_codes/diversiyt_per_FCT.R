setwd("~/OneDrive/University of Helsinki/msc_thesis/viral_species")

Pal6 <- c("#BB3E03", "#E9C46A", "#CA6702", "#0A9396",
          "#005F73", "#80BBA8")

# load in viral physeq obj
physeq <- readRDS("PhySeq_helmi_viral_27.11.23.rds")
# load in viral diversity and richness results
tab <- readRDS("viral_alpha_diversity.rds")

tab_sub <- tab %>% select(diversity_shannon,chao1)
tab_sub$sample <- rownames(tab_sub)

# Load in cluster data
cluster_B4.5 <- read_csv("../bact_comm/HClust_F_CLR_B4.5_18.08.22.csv")
cluster_B4.5$ClusterID <- sub("^", "EC_", cluster_B4.5$ClusterID)
cluster_B7.9 <- read_csv("../bact_comm/HClust_F_CLR_B7.9_18.08.22.csv")
cluster_B7.9$ClusterID <- sub("^", "LC_", cluster_B7.9$ClusterID)
cluster_P <- read_csv("../bact_comm/HClust_F_CLR_P_13.11.23.csv") %>% rename('ClusterID' = "ClusterID_allWGS")
cluster_P$ClusterID <- sub("^", "P_", cluster_P$ClusterID)
cluster_info <- rbind(cluster_B4.5,cluster_B7.9,cluster_P)

plot_data <- left_join(tab_sub,cluster_info, by = c("sample" = "Sample_ID"))
plot_data$Sample_type <- plot_data$sample
plot_data$Sample_type <- gsub(".*_","",plot_data$Sample_type)
plot_data$Sample_type <- gsub('M', 'P', gsub('F', 'P', plot_data$Sample_type)) # comine parental clusters

plot_data_B4 <- plot_data %>% filter(Sample_type == "B4")
plot_data_B5 <- plot_data %>% filter(Sample_type == "B5")
plot_data_B7 <- plot_data %>% filter(Sample_type == "B7")
plot_data_B9 <- plot_data %>% filter(Sample_type == "B9")
plot_data_P <- plot_data %>% filter(Sample_type == "P")
  
comp_1 <- list(c("EC_1","EC_2"),c("EC_2","EC_3"),c("EC_1","EC_3"))
comp_2 <- list(c("LC_1","LC_2"))

# plot viral diversity per FCT
fct_div_plot_per_t <- plot_data_B9 %>% ggplot(aes(x=ClusterID, y=chao1))+
  geom_boxplot()+
  stat_compare_means(comparisons = comp_2, method = "wilcox.test") 
fct_div_plot_per_t

fct_div_plot <- plot_data %>%  ggplot(aes(x=Sample_type, y=diversity_shannon, color=ClusterID))+
  geom_boxplot()+
  theme_bw()+
  scale_color_manual(values=Pal6) +
  theme(legend.position = "none") +
  facet_grid(.~Sample_type, scales = "free_x")
fct_div_plot
ggsave("fct_viral_div_plot.pdf")

fct_rich_plot <- plot_data_B4 %>% ggplot(aes(x=Sample_type, y=chao1, color=ClusterID))+
  geom_boxplot()+
  theme_bw()+
  scale_color_manual(values=Pal6) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = comp_1, method = "wilcox.test") #+
  #facet_grid(.~Sample_type, scales = "free_x")
fct_rich_plot
ggsave("fct_viral_rich_plot.pdf")
