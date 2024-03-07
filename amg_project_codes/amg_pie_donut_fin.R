library(ggplot2)
library(tidyverse)
library(webr)
library(readxl)
library(stringr)
library(gridExtra)

# AMG pie donut multiple plot

setwd("~/OneDrive/University of Helsinki/msc_thesis/amgs")
helmi_vr <- read_xlsx("../viral_selection_2/helmi_catalog_fin.xlsx")
amg_summ <- read_csv("amg_summ_unique.csv") %>% 
  select(gene,scaffold,sheet,header) %>%
  rename("viral_name" = "gene") 

duplicated <- amg_summ$viral_name[duplicated(amg_summ$viral_name)] # select duplicates by scaffold
duplicated <- duplicated %>% unique()
amg_dups <- amg_summ %>% filter(amg_summ$viral_name %in% duplicated) %>% select(-scaffold) # filter to only keep duplicates 

# annotate as multiple
test <- amg_summ
for(i in 1:length(test$viral_name)){
  if(test$viral_name[i] %in% duplicated){
    test$sheet[i] <- "Multiple"
    test$header[i] <- "Multiple"
  }
}

# keep only one multiple annotation per viral name 
amg_summ_mult_derep <- test[!duplicated(test),] 

# annotate remaining C1 as "other"
amg_summ_mult_derep$header <- gsub("C1","Other",amg_summ_mult_derep$header)

# plot 
plot_data <- amg_summ_mult_derep %>% select(-viral_name,-scaffold) %>% group_by(sheet,header) %>% tally()
plot_data$header <- plot_data$header %>% replace_na('Unclassified')

# combine other and photosynthesis energy metabolism 
plot_data[4,3] <- plot_data[4,3] + 16 + 8 + 14 + 1 + 6
plot_data <- plot_data[-c(2,3,5,7,8),]

# combine MISC data
plot_data[4,3] <- plot_data[4,3] + 34 + 28 + 2 + 14 + 19 + 5
plot_data[4,2] <- "Other"
plot_data <- plot_data[-c(5,6,7,9,10,12),]

# combine carbon metabolism
plot_data[13,3] <- plot_data[13,3] + 8 + 107
plot_data[13,2] <- "Other"
plot_data <- plot_data[-c(15,16),]

donut <- PieDonut(plot_data, aes(sheet, header, count=n),labelposition=0,showPieName=FALSE)
donut
ggsave("AMG_donut_15_2_24.PDF", plot = donut)

pie <- PieDonut(plot_data,aes(sheet,count=n),labelpositionThreshold=0.1)
pie
ggsave("AMG_pie_15_2_24.PDF", plot = pie)

prop <- aggregate(plot_data$n, by=list(sheet=plot_data$sheet), FUN=sum)
plot_data <- merge(plot_data,prop,by="sheet")
plot_data$prop <- ((plot_data$n/plot_data$x) * 100)

write.csv(plot_data, file = "amg_plot_data_15_2_24.csv")
