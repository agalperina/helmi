library(ape)
library(dplyr)
library(phyloseq)
library(tidyverse)
library(microbiome)
library(readxl)
library(ggpubr)
library(gridExtra)
library(vegan)
library(plotly)

# 
Pal6 <- c("#BB3E03", "#E9C46A", "#CA6702", "#0A9396",
          "#005F73", "#80BBA8")
Pal6_2<- c("#BB3E03", "#E9C46A", "#CA6702", "#80BBA8", "#0A9396","#005F73")

# Load in data
setwd("~/OneDrive/University of Helsinki/msc_thesis/amgs")
helmi_vr <- read_xlsx("../viral_selection_2/helmi_catalog_fin.xlsx")
amg_summ <- read_csv("amg_summ_unique.csv") %>% 
  select(gene,scaffold,sheet,header) %>%
  rename("viral_name" = "gene") 

amg_summ$scaffold <- gsub("__.*","",amg_summ$scaffold)
duplicated <- amg_summ$scaffold[duplicated(amg_summ$scaffold)] # select duplicates by scaffold
amg_dups <- amg_summ %>% filter(amg_summ$scaffold %in% duplicated) %>% select(-viral_name) # filter to only keep duplicates 

# Need to mark as multiple only the genes that truely have mutiple functions
t <- amg_dups %>% duplicated() 
t2 <- cbind(amg_dups,t)

# determine amgs with mult functions
unique_dups <- t2 %>% filter(t == T) %>% unique() %>% select(scaffold) # duplicated scaffolds appear only one in dataframe 
t2 <- t2 %>% select(-t)
# if there is more than one function type per gene, mark as multiple = 1 if not = 0
for(i in 1:length(unique_dups$scaffold)){
  if(nrow(unique(subset(t2, scaffold == unique_dups$scaffold[i]))) == 1){ # multiple
    unique_dups$mult[i] <- 0
  }else{unique_dups$mult[i] <- 1}
}

# select unique scaffolds b/c only need one to map to abundance tables
amg_summ_u <- amg_summ %>%
  distinct(scaffold, .keep_all=TRUE) 

# change function to multiple if conditions are truw
for(i in 1:length(unique_dups$scaffold)){
  if((unique_dups$scaffold[i] %in% amg_summ_u$scaffold) && (unique_dups[i,2] == 1)){
    for(j in 1:length(amg_summ_u$scaffold)){
      if(amg_summ_u$scaffold[j] == unique_dups$scaffold[i]){
        amg_summ_u[j,3] = "Multiple"
      }
    }
  }
}

# B4
B4_counts <- read_csv2("../viral_species/B4_sel3count_table.csv")
B4_info <- left_join(B4_counts, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
B4_info <- B4_info[,c(1,477,478,2:476)]

B4_info_l <- B4_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B4_info_avg <- B4_info_l %>%
  group_by(Sample) %>% 
  summarize(sum=sum(abd)) 

# function
B4_func_avg <- B4_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

B4_func_avg <- aggregate( .~sheet, data = B4_func_avg, FUN = mean)
B4_func_avg[1,2] <- B4_func_avg[1,2] + B4_func_avg[2,2] # combine carbon utilization (Woodcroft) with carbon utilization
B4_func_avg <- B4_func_avg[-c(2),]
B4_func_avg <- B4_func_avg %>%
  mutate(pro = sum / sum(sum))

#B5 
B5_counts <- read_csv2("../viral_species/B5_sel3count_table.csv")
B5_info <- left_join(B5_counts, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
B5_info <- B5_info[,c(1,477,478,2:476)]

B5_info_l <- B5_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B5_info_avg <- B5_info_l %>%
  group_by(Sample) %>% 
  summarize(sum=sum(abd)) 

# function
B5_func_avg <- B5_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

B5_func_avg <- aggregate( .~sheet, data = B5_func_avg, FUN = mean)
B5_func_avg[1,2] <- B5_func_avg[1,2] + B5_func_avg[2,2] # combine carbon utilization (Woodcroft) with carbon utilization
B5_func_avg <- B5_func_avg[-c(2),]
B5_func_avg <- B5_func_avg %>%
  mutate(pro = sum / sum(sum))

# B7 - use new B7
B7_counts <- read_csv2("B7_sel3count_table.csv")
B7_info <- left_join(B7_counts, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
B7_info <- B7_info[,c(1,477,478,2:476)]

B7_info_l <- B7_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B7_info_avg <- B7_info_l %>%
  group_by(Sample) %>% 
  summarize(sum=sum(abd)) 

# function
B7_func_avg <- B7_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

B7_func_avg <- aggregate( .~sheet, data = B7_func_avg, FUN = mean)
B7_func_avg[1,2] <- B7_func_avg[1,2] + B5_func_avg[2,2] # combine carbon utilization (Woodcroft) with carbon utilization
B7_func_avg <- B7_func_avg[-c(2),]
B7_func_avg <- B7_func_avg %>%
  mutate(pro = sum / sum(sum))

# B9 
B9_counts <- read_csv2("../viral_species/B9_sel3count_table.csv")
B9_info <- left_join(B9_counts, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
B9_info <- B9_info[,c(1,477,478,2:476)]

B9_info_l <- B9_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B9_info_avg <- B9_info_l %>%
  group_by(Sample) %>% 
  summarize(sum=sum(abd)) 

# function
B9_func_avg <- B9_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

B9_func_avg <- aggregate( .~sheet, data = B9_func_avg, FUN = mean)
B9_func_avg[1,2] <- B9_func_avg[1,2] + B9_func_avg[2,2] # combine carbon utilization (Woodcroft) with carbon utilization
B9_func_avg <- B9_func_avg[-c(2),]
B9_func_avg <- B9_func_avg %>%
  mutate(pro = sum / sum(sum))

#M
M_counts <- read_csv2("../viral_species/M_sel3count_table.csv")
M_info <- left_join(M_counts, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
M_info <- M_info[,c(1,306,307,2:305)]

M_info_l <- M_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

M_info_avg <- M_info_l %>%
  group_by(Sample) %>% 
  summarize(sum=sum(abd)) 

# function
M_func_avg <- M_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)

M_func_avg <- aggregate( .~sheet, data = M_func_avg, FUN = mean)
M_func_avg[1,2] <- M_func_avg[1,2] + M_func_avg[2,2] # combine carbon utilization (Woodcroft) with carbon utilization
M_func_avg <- M_func_avg[-c(2),]
M_func_avg <- M_func_avg %>%
  mutate(pro = sum / sum(sum))

#F
F_counts <- read_csv2("../viral_species/F_sel3count_table.csv")
F_info <- left_join(F_counts, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
F_info <- F_info[,c(1,125,126,2:124)]

F_info_l <- F_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

F_info_avg <- F_info_l %>%
  group_by(Sample) %>% 
  summarize(sum=sum(abd)) 

# function
F_func_avg <- F_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  select(-Sample)
  
F_func_avg <- aggregate( .~sheet, data = F_func_avg, FUN = mean)
F_func_avg[1,2] <- F_func_avg[1,2] + F_func_avg[2,2] # combine carbon utilization (Woodcroft) with carbon utilization
F_func_avg <- F_func_avg[-c(2),]
F_func_avg <- F_func_avg %>%
  mutate(pro = sum / sum(sum))

# merge parental data
parental_merge <-rbind(F_func_avg,M_func_avg)[,-3]
parental_merge <-aggregate(sum ~ sheet, data=parental_merge, FUN=sum)
parental_merge <- parental_merge %>%
  mutate(pro = sum / sum(sum))

# combine and plot data

amg_count <- rbind(B4_info_avg,B5_info_avg,B7_info_avg,B9_info_avg,F_info_avg,M_info_avg)
amg_count$Sample <- gsub(".*_","",amg_count$Sample)

amg_count$Sample <- gsub('M', 'P', gsub('F', 'P', amg_count$Sample))

comp <- list(c("B4","B5"),c("B5","B7"), c("B7","B9"),c("B9","P"))
amg_count_plot <- amg_count %>% ggplot(aes(x=Sample, y=sum, color=Sample))+
  geom_boxplot()+
  theme_classic() +
  scale_color_manual(values=Pal6) +
  theme(legend.position = "none") +
  xlab("Timepoints") + ylab("Avg Abundance") +
  stat_compare_means(comparisons = comp, method = "wilcox.test") 
amg_count_plot
#ggsave("amg_count_box_14_2_24.pdf")

amg_func <- rbind(B4_func_avg,B5_func_avg,B7_func_avg,B9_func_avg,parental_merge)
amg_func$sample <- rep(c("B4","B5","B7","B9","P"),times=c(6,6,6,6,6))

amg_func_plot <- amg_func %>%  ggplot(aes(fill=sheet, y=pro, x=sample)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = Pal6_2) +
  theme_classic() +
  xlab("Timepoints") + ylab("Proportion of Avg Abundance") + labs(fill = "Metabolic Function")
amg_func_plot
#ggplotly(amg_func_plot)
#ggsave("amg_func_bar.pdf")

## AMGs in FCTs
cluster_B4.5 <- read_csv("../bact_comm/HClust_F_CLR_B4.5_18.08.22.csv")
cluster_B4.5$ClusterID <- sub("^", "EC_", cluster_B4.5$ClusterID)
cluster_B7.9 <- read_csv("../bact_comm/HClust_F_CLR_B7.9_18.08.22.csv")
cluster_B7.9$ClusterID <- sub("^", "LC_", cluster_B7.9$ClusterID)
cluster_P <- read_csv("../bact_comm/HClust_F_CLR_P_13.11.23.csv") %>% rename('ClusterID' = "ClusterID_allWGS")
cluster_P$ClusterID <- sub("^", "P_", cluster_P$ClusterID)

B4_w_cluster <- left_join(B4_info_avg,cluster_B4.5, by = c("Sample" = "Sample_ID"))
B5_w_cluster <- left_join(B5_info_avg,cluster_B4.5, by = c("Sample" = "Sample_ID"))
B7_w_cluster <- left_join(B7_info_avg,cluster_B7.9, by = c("Sample" = "Sample_ID"))
B9_w_cluster <- left_join(B9_info_avg,cluster_B7.9, by = c("Sample" = "Sample_ID"))
F_w_cluster <- left_join(F_info_avg,cluster_P, by = c("Sample" = "Sample_ID"))
M_w_cluster <- left_join(M_info_avg,cluster_P, by = c("Sample" = "Sample_ID"))

amg_cluster_count <- rbind(B4_w_cluster,B5_w_cluster,B7_w_cluster,B9_w_cluster,F_w_cluster,M_w_cluster)
amg_cluster_count$Sample <- gsub(".*_","",amg_cluster_count$Sample)

amg_cluster_count$Sample <- gsub('M', 'P', gsub('F', 'P', amg_cluster_count$Sample))

amg_cluster_plot <- amg_cluster_count %>% ggplot(aes(x=ClusterID, y=sum, color=ClusterID))+
  geom_boxplot()+
  theme_light()+
  scale_color_manual(values=Pal6) + 
  facet_grid(.~Sample, scales = "free_x") +
  theme(legend.position = "none") +
  xlab("Clusters per Timepoint") + ylab("Avg Abundance")
amg_cluster_plot
#ggsave("amg_cluster_count_box.PDF")

B4_func_cluster <- B4_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd))

B4_func_cluster <- left_join(B4_func_cluster,cluster_B4.5, by = c("Sample" = "Sample_ID"))
B4_func_cluster_avg <- B4_func_cluster %>% 
  group_by(sheet,ClusterID) %>%
  summarise(avg = mean(sum))

# combine carbon utilization (Woodcroft) with carbon utilization
B4_func_cluster_avg[16,3] <- B4_func_cluster_avg[16,3] + B4_func_cluster_avg[19,3] 
B4_func_cluster_avg[17,3] <- B4_func_cluster_avg[17,3] + B4_func_cluster_avg[20,3] 
B4_func_cluster_avg[18,3] <- B4_func_cluster_avg[18,3] + B4_func_cluster_avg[21,3] 

B4_func_cluster_avg <- B4_func_cluster_avg[-c(19,20,21),]
B4_func_cluster_avg <- B4_func_cluster_avg %>%
  group_by(ClusterID) %>%
  mutate(pro = avg / sum(avg))

B5_func_cluster <- B5_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd))

B5_func_cluster <- left_join(B5_func_cluster,cluster_B4.5, by = c("Sample" = "Sample_ID"))
B5_func_cluster_avg <- B5_func_cluster %>% 
  group_by(sheet,ClusterID) %>%
  summarise(avg = mean(sum))

# combine carbon utilization (Woodcroft) with carbon utilization
B5_func_cluster_avg[16,3] <- B5_func_cluster_avg[16,3] + B5_func_cluster_avg[19,3] 
B5_func_cluster_avg[17,3] <- B5_func_cluster_avg[17,3] + B5_func_cluster_avg[20,3] 
B5_func_cluster_avg[18,3] <- B5_func_cluster_avg[18,3] + B5_func_cluster_avg[21,3] 

B5_func_cluster_avg <- B5_func_cluster_avg[-c(19,20,21),]
B5_func_cluster_avg <- B5_func_cluster_avg %>%
  group_by(ClusterID) %>%
  mutate(pro = avg / sum(avg))

B7_func_cluster <- B7_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd))

B7_func_cluster <- left_join(B7_func_cluster,cluster_B7.9, by = c("Sample" = "Sample_ID"))
B7_func_cluster_avg <- B7_func_cluster %>% 
  group_by(sheet,ClusterID) %>%
  summarise(avg = mean(sum))

# combine carbon utilization (Woodcroft) with carbon utilization
B7_func_cluster_avg[11,3] <- B7_func_cluster_avg[11,3] + B7_func_cluster_avg[13,3] 
B7_func_cluster_avg[12,3] <- B7_func_cluster_avg[12,3] + B7_func_cluster_avg[14,3] 

B7_func_cluster_avg <- B7_func_cluster_avg[-c(13,14),]
B7_func_cluster_avg <- B7_func_cluster_avg %>%
  group_by(ClusterID) %>%
  mutate(pro = avg / sum(avg))

B9_func_cluster <- B9_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd))

B9_func_cluster <- left_join(B9_func_cluster,cluster_B7.9, by = c("Sample" = "Sample_ID"))
B9_func_cluster_avg <- B9_func_cluster %>% 
  group_by(sheet,ClusterID) %>%
  summarise(avg = mean(sum))

# combine carbon utilization (Woodcroft) with carbon utilization
B9_func_cluster_avg[11,3] <- B9_func_cluster_avg[11,3] + B9_func_cluster_avg[13,3] 
B9_func_cluster_avg[12,3] <- B9_func_cluster_avg[12,3] + B9_func_cluster_avg[14,3] 

B9_func_cluster_avg <- B9_func_cluster_avg[-c(13,14),]
B9_func_cluster_avg <- B9_func_cluster_avg %>%
  group_by(ClusterID) %>%
  mutate(pro = avg / sum(avg))

M_func_cluster <- M_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd))

M_func_cluster <- left_join(M_func_cluster,cluster_P, by = c("Sample" = "Sample_ID"))
M_func_cluster_avg <- M_func_cluster %>% 
  group_by(sheet,ClusterID) %>%
  summarise(avg = mean(sum))

# combine carbon utilization (Woodcroft) with carbon utilization
M_func_cluster_avg[6,3] <- M_func_cluster_avg[6,3] + M_func_cluster_avg[7,3] 
M_func_cluster_avg <- M_func_cluster_avg[-c(7),]

#M_func_cluster_avg <- M_func_cluster_avg %>%
 # group_by(ClusterID) %>%
 # mutate(pro = avg / sum(avg))

F_func_cluster <- F_info_l %>%
  group_by(sheet, Sample) %>% 
  summarize(sum=sum(abd))

F_func_cluster <- left_join(F_func_cluster,cluster_P, by = c("Sample" = "Sample_ID"))
F_func_cluster_avg <- F_func_cluster %>% 
  group_by(sheet,ClusterID) %>%
  summarise(avg = mean(sum))

# combine carbon utilization (Woodcroft) with carbon utilization
F_func_cluster_avg[6,3] <- F_func_cluster_avg[6,3] + F_func_cluster_avg[7,3] 
F_func_cluster_avg <- F_func_cluster_avg[-c(7),]

#F_func_cluster_avg <- F_func_cluster_avg %>%
 # group_by(ClusterID) %>%
 # mutate(pro = avg / sum(avg))

# merge parental 
parental_merge_2 <- rbind(F_func_cluster_avg, M_func_cluster_avg)
parental_merge_2 <-aggregate(avg ~ sheet, data=parental_merge_2, FUN=mean)
parental_merge_2 <- parental_merge_2 %>%
  mutate(pro = avg / sum(avg)) %>%
  mutate(ClusterID = rep(c("P_1"),times=c(6)))
parental_merge_2 <- parental_merge_2[,c(1,4,2,3)]

amg_cluster_func <- rbind(B4_func_cluster_avg,B5_func_cluster_avg,B7_func_cluster_avg,B9_func_cluster_avg,parental_merge_2)
amg_cluster_func$sample <- rep(c("B4","B5","B7","B9","P"),times=c(18,18,12,12,6))

amg_cluster_func_plot <- amg_cluster_func %>%  ggplot(aes(fill=sheet, y=pro, x=ClusterID)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = Pal6) +
  facet_grid(.~sample, scales = "free_x")+
  theme_bw() 
amg_cluster_func_plot
#ggsave("amg_cluster_func_plot.PDF")

## Functino rel abd 
F_count <- read_csv2("../viral_species/F_sel3count_table.csv")
test <- sweep(F_count[-1], 2, colSums(F_count[,-1]), '/') * 100
virals_name <- F_count[1]
F_count <- cbind(virals_name, test)
F_info <- left_join(F_count, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
F_info <- F_info[,c(1,125,126,2:124)]

F_info_l <- F_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

F_ON <- F_info_l %>% filter(grepl("Organic Nitrogen",sheet)) %>% select(Sample,sheet,abd)
F_carbon <- F_info_l %>% filter(grepl("carbon utilization",sheet)) %>% select(Sample,sheet,abd)


M_count <- read_csv2("../viral_species/M_sel3count_table.csv")
test <- sweep(M_count[-1], 2, colSums(M_count[,-1]), '/') * 100
virals_name <- M_count[1]
M_count <- cbind(virals_name, test)
M_info <- left_join(M_count, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
M_info <- M_info[,c(1,306,307,2:305)]

M_info_l <- M_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

M_ON <- M_info_l %>% filter(grepl("Organic Nitrogen",sheet)) %>% select(Sample,sheet,abd)
M_carbon <- M_info_l %>% filter(grepl("carbon utilization",sheet)) %>% select(Sample,sheet,abd)

B9_count <- read_csv2("../viral_species/B9_sel3count_table.csv")
test <- sweep(B9_count[-1], 2, colSums(B9_count[,-1]), '/') * 100
virals_name <- B9_count[1]
B9_count <- cbind(virals_name, test)
B9_info <- left_join(B9_count, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
B9_info <- B9_info[,c(1,477,478,2:476)]

B9_info_l <- B9_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B9_ON <- B9_info_l %>% filter(grepl("Organic Nitrogen",sheet)) %>% select(Sample,sheet,abd)
B9_carbon <- B9_info_l %>% filter(grepl("carbon utilization",sheet)) %>% select(Sample,sheet,abd)

B7_count <- read_csv2("../viral_species/B7_sel3count_table.csv")
test <- sweep(B7_count[-1], 2, colSums(B7_count[,-1]), '/') * 100
virals_name <- B7_count[1]
B7_count <- cbind(virals_name, test)
B7_info <- left_join(B7_count, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
B7_info <- B7_info[,c(1,477,478,2:476)]

B7_info_l <- B7_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B7_ON <- B7_info_l %>% filter(grepl("Organic Nitrogen",sheet)) %>% select(Sample,sheet,abd)
B7_carbon <- B7_info_l %>% filter(grepl("carbon utilization",sheet)) %>% select(Sample,sheet,abd)

B5_count <- read_csv2("../viral_species/B5_sel3count_table.csv")
test <- sweep(B5_count[-1], 2, colSums(B5_count[,-1]), '/') * 100
virals_name <- B5_count[1]
B5_count <- cbind(virals_name, test)
B5_info <- left_join(B5_count, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
B5_info <- B5_info[,c(1,477,478,2:476)]

B5_info_l <- B5_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B5_ON <- B5_info_l %>% filter(grepl("Organic Nitrogen",sheet)) %>% select(Sample,sheet,abd)
B5_carbon <- B5_info_l %>% filter(grepl("carbon utilization",sheet)) %>% select(Sample,sheet,abd)

B4_count <- read_csv2("../viral_species/B4_sel3count_table.csv")
test <- sweep(B4_count[-1], 2, colSums(B4_count[,-1]), '/') * 100
virals_name <- B4_count[1]
B4_count <- cbind(virals_name, test)
B4_info <- left_join(B4_count, amg_summ_u, by = c("Name" = "scaffold")) %>% drop_na(viral_name)
B4_info <- B4_info[,c(1,477,478,2:476)]

B4_info_l <- B4_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B4_ON <- B4_info_l %>% filter(grepl("Organic Nitrogen",sheet)) %>% select(Sample,sheet,abd)
B4_carbon <- B4_info_l %>% filter(grepl("carbon utilization",sheet)) %>% select(Sample,sheet,abd)

ON_data <- rbind(B4_ON,B5_ON,B7_ON,B9_ON,M_ON,F_ON)
ON_data$Sample <- gsub(".*_","",ON_data$Sample)

ON_data$Sample <- gsub("F","M",ON_data$Sample)

comp <- list(c("B4","M"),c("B5","M"),c("B7","M"),c("B9","M"))
ON_data_plot <- ON_data %>% ggplot(aes(x=Sample, y=abd, color=Sample))+
  geom_boxplot()+
  theme_classic() +
  scale_color_manual(values=Pal6) +
  theme(legend.position = "none") +
  xlab("Timepoints") + ylab("Avg Abundance") +
  stat_compare_means(comparisons = comp, method = "wilcox.test") 
ON_data_plot

carbon_data <- rbind(B4_carbon,B5_carbon,B7_carbon,B9_carbon,M_carbon,F_carbon)
carbon_data$Sample <- gsub(".*_","",carbon_data$Sample)

carbon_data$Sample <- gsub("F","M",carbon_data$Sample)

carbon_data_plot <- carbon_data %>% ggplot(aes(x=Sample, y=abd, color=Sample))+
  geom_boxplot()+
  theme_classic() +
  scale_color_manual(values=Pal6) +
  theme(legend.position = "none") +
  xlab("Timepoints") + ylab("Avg Abundance") +
  stat_compare_means(comparisons = comp, method = "wilcox.test") 
carbon_data_plot



