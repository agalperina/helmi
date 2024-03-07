library(readxl)
library(data.table)  
library(tidyverse)

Pal3 <- c("#CA6702","#0A9396","#005F73")
Pal10 <- c("#9B2226","#AE2012", "#BB3E03", "#CA6702", "#EE9B00", "#E9C46A",
           "#80BBA8", "#0A9396", "#005F73", "#033549")
Pal11 <- c("#9B2226","#AE2012", "#BB3E03", "#CA6702", "#EE9B00", "#E9C46A",
           "#80BBA8", "#0A9396", "#005F73", "#033549","#808080")
Pal12 <- c("#9B2226","#AE2012", "#BB3E03", "#CA6702", "#EE9B00", "#E9C46A",
           "#80BBA8", "#0A9396", "#005F73", "#033549","#808080","#525050")

setwd("~/OneDrive/University of Helsinki/msc_thesis/helmi_iphop/iphop")
# merge iphop results
files <- list.files(pattern = ".csv")
temp <- lapply(files, fread, sep=",")
iphop <- rbindlist(temp)

setwd("~/OneDrive/University of Helsinki/msc_thesis/helmi_iphop")
iphop <- iphop %>% 
  select(-`AAI to closest RaFAH reference`,-`Confidence score`,-`List of methods`) %>%
  separate_wider_delim(`Host genus`, ";", names = c("Host Domain", "Host Phyla", "Host Class","Host Order","Host Family","Host Genus"))

# clean up iphop data
iphop$`Host Domain` <- gsub(".*__","",iphop$`Host Domain`)
iphop$`Host Phyla` <- gsub(".*__","",iphop$`Host Phyla`)
iphop$`Host Class` <- gsub(".*__","",iphop$`Host Class`)
iphop$`Host Order` <- gsub(".*__","",iphop$`Host Order`)
iphop$`Host Family` <- gsub(".*__","",iphop$`Host Family`)
iphop$`Host Genus` <- gsub(".*__","",iphop$`Host Genus`)

iphop <- iphop %>% 
  group_by(Virus) %>%
  summarise(across(everything(), ~paste0(unique(.), collapse = ";"))) %>%
  rename("contig_id" = "Virus")

hbac <- read_excel("../viral_Selection_2/helmi_catalog_fin.xlsx") %>% rename("contig_id" = "viral_name")
amg <- read_csv("../amgs/filt_amg_stat.csv")

mismatch_hbac <- anti_join(hbac,iphop) # missing in iphop but in hbac
mismatch_iphop <- anti_join(iphop, hbac) # missing in hbac but in iphop

#write.csv(mismatch_hbac$contig_id, "hbac_not_iphop.csv")

hbac <- hbac %>% left_join(iphop)
#write.csv(hbac, "helmi_fin_w_iphop.csv")
hbac <- read.csv("helmi_fin_w_iphop.csv")

domain_count <- hbac %>% group_by(`Host.Domain`) %>% tally()

# Determine top 10 bacterial hosts
iphop_fam_count <- hbac %>% group_by(`Host.Family`) %>% tally() %>% top_n(11) #%>% drop_na()
sum(iphop_fam_count$n) # 120558
iphop_fam_count <- iphop_fam_count %>% add_row(`Host.Family` = 'Other',n = (145818 - 120558))

# Keep only top ten bacterial hosts and determine if they have pAMGs
#hbac_top_10 <- hbac %>% 
  #filter(`Host Family` %in% iphop_fam_count$`Host Family`)

hbac_top_10 <- hbac %>% 
  mutate(amg = case_when(
  (hbac$contig_id %in% amg$viral_name) ~ "pAMG present",
  !(hbac$contig_id %in% amg$viral_name) ~ "pAMG absent"
))

# create others category 
hbac_top_10 <- hbac_top_10 %>% 
  mutate(`Host.Family` = case_when(
    !(hbac_top_10$`Host.Family` %in% iphop_fam_count$`Host.Family`) ~ "Other",
    TRUE ~ `Host.Family`
  ))

# plot pAMG proportion
hbac_pro_plot <- hbac_top_10 %>% 
  group_by(`Host.Family`,amg) %>% 
  tally() %>% 
  add_column(total = c(3005,3005,9511,9511,4231,4231,6774,6774,25480,25480,8263,8263,25260,25260,1967,1967,15424,15424,2088,2088,2595,2595,41220,41220)) %>%
  mutate(prop = (n / total) * 100)

hbac_pro_plot$`Host.Family` <- hbac_pro_plot$`Host.Family` #%>% replace_na('Unclassified')

fam_order <- c("Lachnospiraceae","Ruminococcaceae","Bacteroidaceae","Oscillospiraceae","Enterobacteriaceae","Bifidobacteriaceae",
               "Acutalibacteraceae","Veillonellaceae","Streptococcaceae","Rikenellaceae") #,"Other","Unclassified")

ggplot(hbac_pro_plot, aes(x=`Host.Family`, y=n, fill=amg)) +
  geom_bar(stat="identity") + theme_classic() + coord_flip() +
  scale_x_discrete(limits = rev(fam_order)) +
  scale_fill_manual(values= Pal3) +
  labs(x = "Bacterial Family", y = "Nb of Prophages with at least 1 pAMG") 

#ggsave("host_fam_amg.PDF")

# Top ten hosts across sample
hbac_sub_top_10 <- hbac_top_10 %>% select(contig_id,`Host Family`)
hbac_sub_top_10$`Host Family` <- hbac_sub_top_10$`Host Family` %>% replace_na('Unclassified')

cluster_B4.5 <- read_csv("../bact_comm/HClust_F_CLR_B4.5_18.08.22.csv")
cluster_B4.5$ClusterID <- sub("^", "EC_", cluster_B4.5$ClusterID)
cluster_B7.9 <- read_csv("../bact_comm/HClust_F_CLR_B7.9_18.08.22.csv")
cluster_B7.9$ClusterID <- sub("^", "LC_", cluster_B7.9$ClusterID)
cluster_P <- read_csv("../bact_comm/HClust_F_CLR_P_13.11.23.csv") %>% rename('ClusterID' = "ClusterID_allWGS")
cluster_P$ClusterID <- sub("^", "P_", cluster_P$ClusterID)

# B4
B4_counts <- read_csv2("../viral_species/B4_sel3count_table.csv")
B4_info <- left_join(B4_counts, hbac_sub_top_10, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
B4_info <- B4_info[,c(476,1:475)]

B4_info_l <- B4_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B4_fam_avg <- B4_info_l %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_B4.5, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

B4_fam_avg <- aggregate( .~`Host Family` + ClusterID, data = B4_fam_avg, FUN = mean)
B4_fam_avg <- B4_fam_avg %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# B5
B5_counts <- read_csv2("../viral_species/B5_sel3count_table.csv")
B5_info <- left_join(B5_counts, hbac_sub_top_10, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
B5_info <- B5_info[,c(476,1:475)]


B5_info_l <- B5_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B5_fam_avg <- B5_info_l %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_B4.5, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

B5_fam_avg <- aggregate( .~`Host Family`+ ClusterID, data = B5_fam_avg, FUN = mean)
B5_fam_avg <- B5_fam_avg %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# B7
B7_counts <- read_csv2("../viral_species/B7_sel3count_table.csv")
B7_info <- left_join(B7_counts, hbac_sub_top_10, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
B7_info <- B7_info[,c(475,1:474)]

B7_info_l <- B7_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B7_fam_avg <- B7_info_l %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_B7.9, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

B7_fam_avg <- aggregate( .~`Host Family`+ ClusterID, data = B7_fam_avg, FUN = mean)
B7_fam_avg <- B7_fam_avg %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# B9
B9_counts <- read_csv2("../viral_species/B9_sel3count_table.csv")
B9_info <- left_join(B9_counts, hbac_sub_top_10, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
B9_info <- B9_info[,c(476,1:475)]


B9_info_l <- B9_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B9_fam_avg <- B9_info_l %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_B7.9, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

B9_fam_avg <- aggregate( .~`Host Family`+ ClusterID, data = B9_fam_avg, FUN = mean)
B9_fam_avg <- B9_fam_avg %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# M
M_counts <- read_csv2("../viral_species/M_sel3count_table.csv")
M_info <- left_join(M_counts, hbac_sub_top_10, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
M_info <- M_info[,c(305,1:304)]

M_info_l <- M_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

M_fam_avg <- M_info_l %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_P, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

M_fam_avg <- aggregate( .~`Host Family`+ ClusterID, data = M_fam_avg, FUN = mean)

# F
F_counts <- read_csv2("../viral_species/F_sel3count_table.csv")
F_info <- left_join(F_counts, hbac_sub_top_10, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
F_info <- F_info[,c(124,1:123)]


F_info_l <- F_info %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

F_fam_avg <- F_info_l %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_P, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

F_fam_avg <- aggregate( .~`Host Family`+ ClusterID, data = F_fam_avg, FUN = mean)

# merge parental 
parental_merge <- rbind(F_fam_avg, M_fam_avg)
parental_merge <- aggregate(sum ~ `Host Family` + ClusterID, data=parental_merge, FUN=mean)
parental_merge <- parental_merge %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# combine
fam_count <- rbind(B4_fam_avg,B5_fam_avg,B7_fam_avg,B9_fam_avg,parental_merge)
fam_count$sample <- rep(c("B4","B5","B7","B9","P"),times=c(36,36,24,24,12))

fam_count$`Host Family` <- gsub("Unclassified",'zzUnclassified', fam_count$`Host Family`)  
fam_count$`Host Family` <- gsub("Other",'zzOther', fam_count$`Host Family`)  

iphop_fam_plot <- fam_count %>%  ggplot(aes(fill=`Host Family`, y=pro, x=ClusterID)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = Pal12) +
  facet_grid(.~sample, scales = "free_x")+
  theme_bw() +
  xlab("Clusters per Timepoint") + ylab("Avg Abundance")
iphop_fam_plot
ggsave('iphop_fam_plot.pdf')

# Iphop fam plot including unknowns
# Determine top 10 bacterial hosts
iphop_fam_count_w_NA <- hbac %>% group_by(`Host Family`) %>% tally() %>% top_n(11)

# Keep only top ten bacterial hosts and determine if they have pAMGs
hbac_top_10_w_NA <- hbac %>% 
  filter(`Host Family` %in% iphop_fam_count_w_NA$`Host Family`) 

hbac_sub_top_10_NA <- hbac_top_10_w_NA %>% select(contig_id,`Host Family`)

# B4
B4_counts_NA <- read_csv2("../viral_species/B4_sel3count_table.csv")
B4_info_NA <- left_join(B4_counts_NA, hbac_sub_top_10_NA, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
B4_info_NA <- B4_info_NA[,c(476,1:475)]

B4_info_l_NA <- B4_info_NA %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B4_fam_avg_NA <- B4_info_l_NA %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_B4.5, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample) 

B4_fam_avg_NA$`Host Family` <- B4_fam_avg_NA$`Host Family` %>% replace_na('Unclassified')

B4_fam_avg_NA <- aggregate( .~`Host Family` + ClusterID, data = B4_fam_avg_NA, FUN = mean)
B4_fam_avg_NA <- B4_fam_avg_NA %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# B5
B5_counts_NA <- read_csv2("../viral_species/B5_sel3count_table.csv")
B5_info_NA <- left_join(B5_counts_NA, hbac_sub_top_10_NA, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
B5_info_NA <- B5_info_NA[,c(476,1:475)]

B5_info_l_NA <- B5_info_NA %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B5_fam_avg_NA <- B5_info_l_NA %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_B4.5, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

B5_fam_avg_NA$`Host Family` <- B5_fam_avg_NA$`Host Family` %>% replace_na('Unclassified')

B5_fam_avg_NA <- aggregate( .~`Host Family`+ ClusterID, data = B5_fam_avg_NA, FUN = mean)
B5_fam_avg_NA <- B5_fam_avg_NA %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# B7
B7_counts_NA <- read_csv2("../viral_species/B7_sel3count_table.csv")
B7_info_NA <- left_join(B7_counts_NA, hbac_sub_top_10_NA, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
B7_info_NA <- B7_info_NA[,c(475,1:474)]

B7_info_l_NA <- B7_info_NA%>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B7_fam_avg_NA <- B7_info_l_NA %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_B7.9, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

B7_fam_avg_NA$`Host Family` <- B7_fam_avg_NA$`Host Family` %>% replace_na('Unclassified')

B7_fam_avg_NA <- aggregate( .~`Host Family`+ ClusterID, data = B7_fam_avg_NA, FUN = mean)
B7_fam_avg_NA <- B7_fam_avg_NA %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# B9
B9_counts_NA <- read_csv2("../viral_species/B9_sel3count_table.csv")
B9_info_NA <- left_join(B9_counts_NA, hbac_sub_top_10_NA, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
B9_info_NA <- B9_info_NA[,c(476,1:475)]


B9_info_l_NA <- B9_info_NA %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

B9_fam_avg_NA <- B9_info_l_NA %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_B7.9, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

B9_fam_avg_NA$`Host Family` <- B9_fam_avg_NA$`Host Family` %>% replace_na('Unclassified')

B9_fam_avg_NA <- aggregate( .~`Host Family`+ ClusterID, data = B9_fam_avg_NA, FUN = mean)
B9_fam_avg_NA <- B9_fam_avg_NA %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# M
M_counts_NA <- read_csv2("../viral_species/M_sel3count_table.csv")
M_info_NA <- left_join(M_counts_NA, hbac_sub_top_10_NA, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
M_info_NA <- M_info_NA[,c(305,1:304)]

M_info_l_NA <- M_info_NA %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

M_fam_avg_NA <- M_info_l_NA %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_P, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

M_fam_avg_NA$`Host Family` <- M_fam_avg_NA$`Host Family` %>% replace_na('Unclassified')

M_fam_avg_NA <- aggregate( .~`Host Family`+ ClusterID, data = M_fam_avg_NA, FUN = mean)

# F
F_counts_NA <- read_csv2("../viral_species/F_sel3count_table.csv")
F_info_NA <- left_join(F_counts_NA, hbac_sub_top_10_NA, by = c("Name" = "contig_id")) %>% column_to_rownames(var="Name")
F_info_NA <- F_info_NA[,c(124,1:123)]


F_info_l_NA <- F_info_NA %>%
  pivot_longer(
    cols = starts_with("Perhe"),
    names_to = "Sample",
    values_to = "abd"
  ) 

F_fam_avg_NA <- F_info_l_NA %>%
  group_by(`Host Family`, Sample) %>% 
  summarize(sum=sum(abd)) %>%
  left_join(cluster_P, by = c("Sample" = "Sample_ID")) %>%
  select(-Sample)

F_fam_avg_NA$`Host Family` <- F_fam_avg_NA$`Host Family` %>% replace_na('Unclassified')

F_fam_avg_NA <- aggregate( .~`Host Family`+ ClusterID, data = F_fam_avg_NA, FUN = mean)

# merge parental 
parental_merge_NA <- rbind(F_fam_avg_NA, M_fam_avg_NA)
parental_merge_NA <- aggregate(sum ~ `Host Family` + ClusterID, data=parental_merge_NA, FUN=mean)
parental_merge_NA <- parental_merge_NA %>%
  group_by(ClusterID) %>%
  mutate(pro = sum / sum(sum))

# combine
fam_count_NA <- rbind(B4_fam_avg_NA,B5_fam_avg_NA,B7_fam_avg_NA,B9_fam_avg_NA,parental_merge_NA)
fam_count_NA$sample <- rep(c("B4","B5","B7","B9","P"),times=c(33,33,22,22,11))

fam_count_NA$`Host Family` <- gsub("Unclassified","zzUnclassified",fam_count_NA$`Host Family`)

iphop_fam_plot_NA <- fam_count_NA %>%  ggplot(aes(fill=`Host Family`, y=pro, x=ClusterID)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = Pal11) +
  facet_grid(.~sample, scales = "free_x")+
  theme_bw() +
  xlab("Clusters per Timepoint") + ylab("Avg Abundance")
iphop_fam_plot_NA
ggsave("iphop_family_plot.PDF")
plotly::ggplotly(iphop_fam_plot_NA)
