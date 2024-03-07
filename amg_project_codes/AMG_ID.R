library(ggplot2)
library(tidyverse)
library(webr)
library(readxl)
library(stringr)
library(gridExtra)

setwd("~/OneDrive/University of Helsinki/msc_thesis/amgs")
Pal3 <- c("#CA6702","#0A9396","#005F73")

# load in data 
amg_summ_a <- read_tsv("dramv_out/xaa_amg_summary.tsv")
amg_summ_b <- read_tsv("dramv_out/xab_amg_summary.tsv")
amg_summ_c <- read_tsv("dramv_out/xac_amg_summary.tsv")
amg_summ_d <- read_tsv("dramv_out/xad_amg_summary.tsv")
amg_summ_e <- read_tsv("dramv_out/xae_amg_summary.tsv")
amg_summ_f <- read_tsv("dramv_out/xaf_amg_summary.tsv")
amg_summ_g <- read_tsv("dramv_out/xag_amg_summary.tsv")
amg_summ_h <- read_tsv("dramv_out/xah_amg_summary.tsv")
amg_summ_i <- read_tsv("dramv_out/xai_amg_summary.tsv")
amg_summ_j <- read_tsv("dramv_out/xaj_amg_summary.tsv")
amg_summ_k <- read_tsv("dramv_out/xak_amg_summary.tsv")
amg_summ_l <- read_tsv("dramv_out/xal_amg_summary.tsv")
amg_summ_m <- read_tsv("dramv_out/xam_amg_summary.tsv")
amg_summ_n <- read_tsv("dramv_out/xan_amg_summary.tsv")
amg_summ_o <- read_tsv("dramv_out/xao_amg_summary.tsv")

amg_stat_a <- read_tsv("dramv_out/xaa_vMAG_stats.tsv")
amg_stat_b <- read_tsv("dramv_out/xab_vMAG_stats.tsv")
amg_stat_c <- read_tsv("dramv_out/xac_vMAG_stats.tsv")
amg_stat_d <- read_tsv("dramv_out/xad_vMAG_stats.tsv")
amg_stat_e <- read_tsv("dramv_out/xae_vMAG_stats.tsv")
amg_stat_f <- read_tsv("dramv_out/xaf_vMAG_stats.tsv")
amg_stat_g <- read_tsv("dramv_out/xag_vMAG_stats.tsv")
amg_stat_h <- read_tsv("dramv_out/xah_vMAG_stats.tsv")
amg_stat_i <- read_tsv("dramv_out/xai_vMAG_stats.tsv")
amg_stat_j <- read_tsv("dramv_out/xaj_vMAG_stats.tsv")
amg_stat_k <- read_tsv("dramv_out/xak_vMAG_stats.tsv")
amg_stat_l <- read_tsv("dramv_out/xal_vMAG_stats.tsv")
amg_stat_m <- read_tsv("dramv_out/xam_vMAG_stats.tsv")
amg_stat_n <- read_tsv("dramv_out/xan_vMAG_stats.tsv")
amg_stat_o <- read_tsv("dramv_out/xao_vMAG_stats.tsv")

# merge 
amg_sum <- rbind(amg_summ_a,amg_summ_b,amg_summ_c,amg_summ_d,amg_summ_e,amg_summ_f,amg_summ_g,amg_summ_h,amg_summ_i,amg_summ_j,
              amg_summ_k,amg_summ_l,amg_summ_m,amg_summ_n,amg_summ_o)

amg_stat <- rbind(amg_stat_a,amg_stat_b,amg_stat_c,amg_stat_d,amg_stat_e,amg_stat_f,amg_stat_g,amg_stat_h,amg_stat_i,amg_stat_j,
                  amg_stat_k,amg_stat_l,amg_stat_m,amg_stat_n,amg_stat_o)

rm(amg_summ_a,amg_summ_b,amg_summ_c,amg_summ_d,amg_summ_e,amg_summ_f,amg_summ_g,amg_summ_h,amg_summ_i,amg_summ_j,
   amg_summ_k,amg_summ_l,amg_summ_m,amg_summ_n,amg_summ_o,amg_stat_a,amg_stat_b,amg_stat_c,amg_stat_d,amg_stat_e,amg_stat_f,amg_stat_g,amg_stat_h,amg_stat_i,amg_stat_j,
   amg_stat_k,amg_stat_l,amg_stat_m,amg_stat_n,amg_stat_o)

# Filter out the uncertain AMGs
# Filter out AMGs with flags

non_amg <- amg_sum %>% filter(potential_amg == TRUE)
amg_E <- non_amg %>% filter((grepl("F",amg_flags)))

AMG_summ_filt <- amg_sum %>% 
  filter(potential_amg == TRUE) %>%
  filter(!(grepl("F",amg_flags)))

AMG_summ_filt <- AMG_summ_filt %>% select(gene,gene_id,scaffold,sheet,header) %>% unique()

##################### This plotting data is outdated, new pie-donut plot is in amg_pie_donut_fin.R file ##################### 
# plot does not include "multiple" category 
# Creating plotting data
total <- read_csv("amg_total_anno.csv")[,-1]
total_plot <- read_csv("amg_plot_anno.csv")[,-1]

total_plot_prop <- c((367/532),(165/532),(217/6536),(6099/6536),(220/6536),(7675/9030),(1355/9030),(1592/1592),(506/1502),(468/1502),(234/1502),(294/1502))

donut <- PieDonut(total_plot, aes(sheet, header, count=n),labelposition=0,showPieName=FALSE)
donut
ggsave("AMG_donut.PDF", plot = donut)

pie <- PieDonut(total_plot,aes(sheet,count=n),labelpositionThreshold=0.1)
pie
ggsave("AMG_pie.PDF", plot = pie)

##############################################################################################################################   

# load in viral catalog 
helmi_vr <- read_csv("../helmi_iphop/helmi_fin_w_iphop.csv")

# How many seq carry at least one AMG
amg_0 <- amg_stat %>% filter(`potential AMG count` == 0)
AMG_1 <- amg_stat %>% filter(`potential AMG count` > 0) %>% 
  rename("contig_id" = "...1") 

# filter out from stat the seq that do not appear in summary 
AMG_1 <- AMG_1[(AMG_1$contig_id %in% AMG_summ_filt$scaffold),]
# 122236  have at least 1 AMG with some carrying up to 43 AMGs

AMG_1$contig_id <- str_remove(AMG_1$contig_id, "__.*")

# check that all AMGs are viral
helmi_sub <- helmi_vr %>%
  select(contig_id,genomad_Class, phagcn_Family, lifestyle, `Host Genus`)

AMG_1_join <- left_join(AMG_1, helmi_sub, by = "contig_id")
levels(factor(AMG_1_join$genomad_Class)) # yes only phages

write_csv(AMG_1_join,"filt_amg_stat.csv")
write_csv(AMG_summ_filt,"filt_amg_summ.csv")


# Calculate the proportion of viral sequences with AMGs 
# according to the viral families (predicted by Genomad/PhaGCN) 

AMG_1_Class <- AMG_1_join %>% group_by(genomad_Class) %>% tally() 
AMG_1_Family <- AMG_1_join %>% group_by(genomad_Class,phagcn_Family) %>% tally() 
AMG_1_lifestyle <- AMG_1_join %>% group_by(lifestyle) %>% tally() 
AMG_1_lifestyle_2 <- AMG_1_join %>% group_by(phagcn_Family,lifestyle) %>% tally() 
#AMG_1_lifestyle_3 <- AMG_1_join %>% group_by(`Host Genus`,lifestyle) %>% tally() %>% filter(`Host Genus` == "Lachnospira" | `Host Genus` == "Bacteroides" | `Host Genus` == "Bifidobacterium")

class <- helmi_vr %>% group_by(genomad_Class) %>% tally()
fam <- helmi_vr %>% group_by(genomad_Class,phagcn_Family) %>% tally()

AMG_1_Family <- cbind(AMG_1_Family, fam[1:21,3])
colnames(AMG_1_Family)[3] <- "n"
colnames(AMG_1_Family)[4] <- "total"
AMG_1_Family$diff <- (AMG_1_Family$total - AMG_1_Family$n)
AMG_1_Family$prop <- ((AMG_1_Family$n / AMG_1_Family$total) * 100)

fam_plot <- data_frame(Family = c("Ackermannviridae","Ackermannviridae","Autographiviridae","Autographiviridae","Casjensviridae","Casjensviridae",
                               "Chaseviridae","Chaseviridae","Demerecviridae","Demerecviridae","Drexlerviridae","Drexlerviridae","Guelinviridae","Guelinviridae",
                               "Herelleviridae","Herelleviridae","Kyanoviridae","Kyanoviridae","Mesyanzhinovviridae","Mesyanzhinovviridae",
                               "Orlajensenviridae","Orlajensenviridae","Peduoviridae","Peduoviridae","Rountreeviridae","Rountreeviridae",
                               "Salasmaviridae","Salasmaviridae","Schitoviridae","Schitoviridae","Straboviridae","Straboviridae",
                               "Vilmaviridae","Vilmaviridae","Zierdtviridae","Zierdtviridae","Zobellviridae","Zobellviridae"),
                       Proportion = c ("AMG present","AMG absent","AMG present","AMG absent","AMG present","AMG absent",
                                       "AMG present","AMG absent","AMG present","AMG absent","AMG present","AMG absent","AMG present","AMG absent",
                                       "AMG present","AMG absent","AMG present","AMG absent","AMG present","AMG absent",
                                       "AMG present","AMG absent","AMG present","AMG absent","AMG present","AMG absent",
                                       "AMG present","AMG absent","AMG present","AMG absent","AMG present","AMG absent",
                                       "AMG present","AMG absent","AMG present","AMG absent","AMG present","AMG absent"),
                       n = c(603,3417,109,273,2172,8107,311,995,146,253,524,3825,127,598,924,2593,63,477,151,2272,32,744,1750,13759,3,86,444,3206,
                             56,310,2851,11015,20,527,104,629,3,18))

fam_order <- c("Peduoviridae","Straboviridae","Casjensviridae","Drexlerviridae","Ackermannviridae","Salasmaviridae",
               "Herelleviridae","Mesyanzhinovviridae","Chaseviridae","Orlajensenviridae","Zierdtviridae","Guelinviridae","Vilmaviridae",
               "Kyanoviridae","Demerecviridae","Autographiviridae","Schitoviridae","Rountreeviridae","Zobellviridae")

ggplot(fam_plot, aes(x=Family, y=n, fill=Proportion)) +
  geom_bar(stat="identity") + theme_classic() + coord_flip() +
  scale_x_discrete(limits = rev(fam_order)) +
  scale_fill_manual(values= Pal3) +
  labs(x = "Number AMGs", y = "Family") +
  ggtitle("Proportion AMGs in Caudoviricetes Families") +
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),legend.title = element_text(size = 14), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14))

test <- aggregate(fam_plot$n , by=list(Family=contig_number$Family), FUN=sum)
contig_number_2 <-  merge(contig_number,test,by="Family")
#ggsave("amg_fam_dist.pdf")

################################ Creating plotting data ################################

# Summarize the main functions categories found as AMG
total <- AMG_summ_filt %>% group_by(sheet,header) %>% tally() 
total$header <- total$header %>% replace_na('Unclassified')

pdf("amg_metabolic_table.pdf", height=8, width=10)
grid.table(total)
dev.off()

total_plot <- total

# Cleaning up total data
total_plot[25,3] <- total_plot[25,3] + 272
total_plot <- total_plot[-c(27),]

# keep top 3
total_plot[2,3] <- 165
total_plot[2,2] <- "Other"
total_plot <- total_plot[-c(3:8),]

total_plot[3,3] <- 25 + 71 + 97 + 19 + 5
total_plot[3,2] <- "Other"
total_plot <- total_plot[-c(4:6,8,9,11),]

total_plot[11,3] <- 40 + 186 + 8
total_plot[11,2] <- "Other"
total_plot <- total_plot[-c(12,14),]

write.csv(total, file = "amg_total_anno.csv")
write.csv(total_plot, file = "amg_plot_anno.csv")


