library(readxl)
library(data.table)  
library(tidyverse)

# file where we merge avrc and iphop

Pal2 <- c("#CA6702","#005F73")
Pal3 <- c("#CA6702","#0A9396","#005F73")
Pal5 <- c("#CA6702", "#E9C46A","#0A9396","#005F73","#80BBA8")
Pal6 <- c("#CA6702","#BB3E03", "#E9C46A","#0A9396","#005F73","#80BBA8")
Pal11 <- c("#9B2226","#AE2012", "#BB3E03", "#CA6702", "#EE9B00", "#E9C46A",
           "#80BBA8", "#0A9396", "#005F73", "#033549","#808080")
Pal20 <- c("#FFBA08","#FAA307","#F48C06","#E85D04","#DC2F02","#D00000","#9D0208","#6A040F","#370617","#03071E","#012A4A","#013A63", "#01497C", "#014F86", "#2A6F97", "#2C7DA0",
           "#468FAF", "#61A5C2", "#89C2D9", "#808080")
Pal21 <- c("#FFBA08","#FAA307","#F48C06","#E85D04","#DC2F02","#D00000","#9D0208","#6A040F","#370617","#03071E","#012A4A","#013A63", "#01497C", "#014F86", "#2A6F97", "#2C7DA0",
           "#468FAF", "#61A5C2", "#89C2D9", "#caf0f8", "#808080")
setwd("~/OneDrive/University of Helsinki/msc_thesis/avrc")
avrc <- read_excel("unified_catalog_final.xlsx")

# AVrC characteristics 
# quality
avrc_quality <- avrc
avrc_quality$study <- avrc_quality$contig_id
avrc_quality$study <- gsub("GutCatV1_","",avrc_quality$study)
avrc_quality$study <- gsub("_.*","",avrc_quality$study)

avrc_quality_p <- avrc_quality %>% group_by(study,checkv_quality) %>% tally()
ggplot(avrc_quality_p, aes(fill=checkv_quality, y=n, x=study)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values= Pal3) +
  theme_classic() +
  labs(x = "Study", y = "Value")
#ggsave("avrc_quality.pdf")

# lifestyle
avrc_lifestyle <- avrc_quality %>% group_by(study,lifestyle) %>% tally()
ggplot(avrc_lifestyle, aes(fill=lifestyle, y=n, x=study)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values= Pal5) +
  theme_classic() +
  labs(x = "Study", y = "Value")
#ggsave("avrc_lifestyle.pdf")

# class distribution
avrc_class <- avrc_quality %>% group_by(genomad_Class) %>% tally()
ggplot(avrc_class, aes(x="", y=n, fill=genomad_Class)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values= Pal20) +
  theme_classic()
ggsave("avrc_class.pdf") 

# merging iphop together
# merge to gether iphop results
setwd("~/OneDrive/University of Helsinki/msc_thesis/avrc/iphop")
files <- list.files(pattern = ".csv")
temp <- lapply(files, fread, sep=",")
iphop <- rbindlist(temp)

# split tax into many columns
iphop <- iphop %>%
separate_wider_delim(`Host genus`, ";", names = c("Host Domain", "Host Phyla", "Host Class","Host Order","Host Family","Host Genus"))
iphop$`Host Domain` <- gsub(".*__","",iphop$`Host Domain`)
iphop$`Host Phyla` <- gsub(".*__","",iphop$`Host Phyla`)
iphop$`Host Class` <- gsub(".*__","",iphop$`Host Class`)
iphop$`Host Order` <- gsub(".*__","",iphop$`Host Order`)
iphop$`Host Family` <- gsub(".*__","",iphop$`Host Family`)
iphop$`Host Genus` <- gsub(".*__","",iphop$`Host Genus`)

#write.csv(iphop, "iphop_merged.csv")
iphop <- read_csv("iphop_merged.csv")
iphop_p <- iphop
iphop_p$study <- iphop$Virus
iphop_p$study <- gsub("GutCatV1_","",iphop_p$study)
iphop_p$study <- gsub("_.*","",iphop_p$study)

iphop_p_2 <- iphop_p %>% group_by(`Host Class`) %>% tally() %>% top_n(10) %>% add_row(`Host Class` = "zzOther", n = 199317 - 195957)
sum(iphop_p_2$n)

# plot iphop results
ggplot(iphop_p_2, aes(x="", y=n, fill=`Host Class`)) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values= Pal11) +
  theme_classic()
ggsave("avrc_iphop_host_13.2.24.pdf")

# 
# Only appear in one catalog sequences
avrc_6 <- avrc %>%
  select(GVDv1:Koonin)

to_integer <- function(x){
  as.integer(as.logical(x))
}
avrc_6 <- avrc_6 %>% 
  mutate(across(GVDv1:Koonin, to_integer),sum = rowSums(.)) 

avrc_6 <- cbind(avrc$contig_id,avrc_6)
avrc_unique <- avrc_6 %>% filter(sum == 1) %>% rename("contig_id" = "avrc$contig_id")

avrc_unique_plot <- avrc_quality %>% left_join(avrc_unique) %>% group_by(study,sum) %>% tally() %>% replace(is.na(.), 0)
as.character(avrc_unique_plot$sum)
avrc_unique_plot$sum <- gsub("1","Unique",avrc_unique_plot$sum)
avrc_unique_plot$sum <- gsub("0","Non-unique",avrc_unique_plot$sum)

ggplot(avrc_unique_plot, aes(fill=sum, y=n, x=study)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values= Pal2, name = "Overlap") +
  theme_classic() +
  labs(x = "Study", y = "Value") 
#ggsave("avrc_overlap.pdf")

# merging avrc with iphop data 
iphop <- read_csv("iphop_merged.csv")[,-1]
iphop <- iphop %>% select(-`AAI to closest RaFAH reference`,-`Confidence score`,-`List of methods`)

iphop <- iphop %>% 
  group_by(Virus) %>%
  summarise(across(everything(), ~paste0(unique(.), collapse = ";"))) %>%
  rename("contig_id" = "Virus")

avrc <- avrc %>% left_join(iphop)

test <- avrc[!avrc$contig_id %in% iphop$contig_id,]

# antijoin what is not in avrc 
mismatch_avrc <- anti_join(avrc,iphop) # missing in iphop but in avrc
mismatch_iphop <- anti_join(iphop, avrc) # missing in avrc but in iphop

write.csv(mismatch_avrc, "avrc_not_iphop.csv")
write.csv(mismatch_iphop, "iphop_not_avrc.csv")

# non iphop annotated in avrc

# most common domain 
avrc_domain <- avrc %>% group_by(`Host Domain`) %>% tally()
avrc_class_bact <- avrc %>% filter(`Host Domain` == "Bacteria") %>% group_by(`Host Class`) %>% tally()  
avrc_tot_classes <-  avrc_class_bact %>% filter(!str_detect(`Host Class`, ';')) # + Desulfurobacteriia + Polyangia + Rhodothermia + Phycisphaerae  

write.csv(avrc, "avrc_w_iphop.csv")
