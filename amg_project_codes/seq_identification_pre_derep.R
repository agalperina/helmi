library(ggplot2)
library(tidyverse)
library(dplyr)
library(grid)
library(gridExtra)
library(plotly)
library(xlsx)
library(ggpubr)
library(stats)
library(corrplot) 

setwd("~/OneDrive/University of Helsinki/msc_thesis/viral_seq_id_pre_derep")

# Colors 
Pal3 <- c("#CA6702","#0A9396","#005F73")
Pal5 <- c("#AE2012", "#E9C46A", "#CA6702","#80BBA8", "#005F73")
Pal5_2 <- c("#CA6702","#AE2012","#80BBA8", "#E9C46A","#005F73")


Pal6 <- c("#BB3E03", "#E9C46A", "#CA6702", "#0A9396",
          "#005F73", "#80BBA8")

Pal10 <- c("#9B2226","#AE2012", "#BB3E03", "#CA6702", "#EE9B00", "#E9C46A",
           "#80BBA8", "#0A9396", "#005F73", "#033549")

##

# Load in viral data
#phage_table <- read_csv("VP2_all_selection2.csv")
phage_table <- read_csv("helmi_for_plotting.csv")[,c(-1)]


# Split data by time point
sample_type <- phage_table %>% group_by(Sample) 
sample_type$Sample <- gsub('M', 'P',
                           gsub('F', 'P', sample_type$Sample))
sample_type_2 <- group_split(sample_type) %>% setNames(c("B4", "B5", "B7", "B9", "P"))

# Determine which technology returned information
phage_B4 <- sample_type_2$B4 %>% mutate(Tools = ifelse((DVF & VS2),"Both",
                                                     ifelse(DVF, "DVF", "VS2")))
phage_B5 <- sample_type_2$B5 %>% mutate(Tools = ifelse((DVF & VS2),"Both",
                                                     ifelse(DVF, "DVF", "VS2")))
phage_B7 <- sample_type_2$B7 %>% mutate(Tools = ifelse((DVF & VS2),"Both",
                                                     ifelse(DVF, "DVF", "VS2")))
phage_B9 <- sample_type_2$B9 %>% mutate(Tools = ifelse((DVF & VS2),"Both",
                                                     ifelse(DVF, "DVF", "VS2")))
phage_P <- sample_type_2$P %>% mutate(Tools = ifelse((DVF & VS2),"Both",
                                                   ifelse(DVF, "DVF", "VS2")))


# Set up data for plotting
tool_B4 <- phage_B4 %>% count(Tools) %>%  mutate(time = rep(c("12 weeks"),each=3))
tool_B5 <- phage_B5 %>% count(Tools) %>%  mutate(time = rep(c("6 months"),each=3))
tool_B7 <- phage_B7 %>% count(Tools) %>%  mutate(time = rep(c("12 months"),each=3))
tool_B9 <- phage_B9 %>% count(Tools) %>%  mutate(time = rep(c("24 months"),each=3))
tool_P <- phage_P %>% count(Tools) %>%  mutate(time = rep(c("Parental"),each=3))

# Plotting data
tool_graphing_data <- bind_rows(tool_B4, tool_B5, tool_B7, tool_B9,tool_P)
# To correct that r puts x axis in alphabetial order 
level_order <- c('12 weeks', '6 months', '12 months', '24 months', 'Parental') 

# Plot
tool_graph <- ggplot(tool_graphing_data, aes(fill=Tools, y=n, x=time)) + 
  geom_bar(position="fill", stat="identity") + 
  labs(x = "Age at Sample Collection", y = "Proportion") +
  scale_fill_manual(values= Pal3) +
  theme_classic() +
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 18),legend.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.text = element_text(size = 18)) + 
  scale_x_discrete(limits = level_order) 
tool_graph
ggsave("tool_graph.pdf")

## CheckV quality 

checkV_B4 <- phage_B4 %>% count(checkv_quality) %>%  mutate(Time = rep("12 weeks",5))
checkV_B5 <- phage_B5 %>% count(checkv_quality) %>%  mutate(Time = rep("6 months",5))
checkV_B7 <- phage_B7 %>% count(checkv_quality) %>%  mutate(Time = rep("12 months",5))
checkV_B9 <- phage_B9 %>% count(checkv_quality) %>%  mutate(Time = rep("24 months",5))
checkV_P <- phage_P %>% count(checkv_quality) %>%  mutate(Time = rep("Parental",5))

checkV_graphing_data <- bind_rows(checkV_B4, checkV_B5, checkV_B7, checkV_B9, checkV_P)

# pivot wider to create table for fisher test
checkV_graphing_data_w <- checkV_graphing_data %>%
  pivot_wider(names_from = Time, values_from = n) %>% 
  column_to_rownames(var="checkv_quality")

#write_csv(checkV_graphing_data_w,"checkV_graphing_data_w.csv")
#write_csv(checkV_graphing_data,"checkV_graphing_data.csv")

chisq <- chisq.test(checkV_graphing_data_w) 

corrplot(chisq$residuals, is.cor = FALSE)

checkV_graphing_data$checkv_quality <- factor(checkV_graphing_data$checkv_quality , levels=c('Complete','High-quality','Medium-quality','Low-quality','Not-determined'))
checkV_graph <- ggplot(checkV_graphing_data, aes(fill=checkv_quality, y=n, x=Time)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Age at Sample Collection", y = "Proportion") +
  scale_fill_manual(values=Pal6) + 
  theme_classic() +
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 18),legend.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  guides(fill=guide_legend(title="Quality")) +
  scale_x_discrete(limits = level_order) 
checkV_graph
ggsave("checkV_13.2.24.pdf")

## Boxplot showing the length of viral sequences per type of sample 

contig_leng_B4 <- phage_B4 %>% select(Sample,contig_length)
contig_leng_B5 <- phage_B5 %>% select(Sample,contig_length)
contig_leng_B7 <- phage_B7 %>% select(Sample,contig_length)
contig_leng_B9 <- phage_B9 %>% select(Sample,contig_length)
contig_leng_P <- phage_P %>% select(Sample,contig_length)

contig_graphing_data <- bind_rows(contig_leng_B4, contig_leng_B5, contig_leng_B7, contig_leng_B9, contig_leng_P)
# Have to log transform contig lengths for graph because otherwise they are too long and skew the graph to not
# being able to see the boxplot
contig_graphing_data$Sample <- gsub("B4", "12 weeks", contig_graphing_data$Sample)
contig_graphing_data$Sample <- gsub("B5", "6 months", contig_graphing_data$Sample)
contig_graphing_data$Sample <- gsub("B7", "12 months", contig_graphing_data$Sample)
contig_graphing_data$Sample <- gsub("B9", "24 months", contig_graphing_data$Sample)
contig_graphing_data$Sample <- gsub("P", "Parental", contig_graphing_data$Sample)
contig_graphing_data$contig_length <- log(contig_graphing_data$contig_length)

comp <- list(c("12 weeks","6 months"),c("6 months","12 months"),c("12 months","24 months"),c("24 months","Parental"))
contig_graph <- ggplot(contig_graphing_data, aes(x=Sample, y=contig_length, color=Sample)) +
  geom_boxplot(outlier.size = 0.01) +
  labs(x = "Age at Sample Collection", y = "Log Contig Length") +
  scale_colour_manual(values=Pal5_2)+
  theme_classic() +
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 18),legend.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.text = element_text(size = 18), legend.position = "none") +
  stat_compare_means(comparisons = comp, method = "wilcox.test") +
  scale_x_discrete(limits = level_order) 
contig_graph
ggsave("contig_graph_14_2_24.pdf")

# Number of contigs per sample
contig_number <- phage_table %>% count(Family,Sample)
contig_number$Sample <- gsub("B4", "12 weeks", contig_number$Sample)
contig_number$Sample <- gsub("B5", "6 months", contig_number$Sample)
contig_number$Sample <- gsub("B7", "12 months", contig_number$Sample)
contig_number$Sample <- gsub("B9", "24 months", contig_number$Sample)
contig_number$Sample <- gsub('M', 'Parental',
                             gsub('F', 'Parental', contig_number$Sample))

test <- aggregate(contig_number$n, by=list(Family=contig_number$Family), FUN=sum)
contig_number_2 <-  merge(contig_number,test,by="Family")
contig_number_2$mill <- ((contig_number_2$n / contig_number_2$x) * 1000000)

contig_number_graph <- ggplot(contig_number_2, aes(x=Sample, y=mill, color=Sample)) +
  geom_boxplot(outlier.size = 0.01) +
  labs(x = "Age at Sample Collection", y = "Number of Contigs") +
  scale_colour_manual(values=Pal5_2)+
  theme_classic() +
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 18),legend.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.text = element_text(size = 18), legend.position = "none") +
  stat_compare_means(comparisons = comp, method = "wilcox.test") +
  #guides(fill=guide_legend(title="Sample Type"))  +
  scale_x_discrete(limits = level_order) 
contig_number_graph
ggsave("contig_number_graph_14_2_24.pdf")





