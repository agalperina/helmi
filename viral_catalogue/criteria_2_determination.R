library(ggplot2)
library(tidyverse)
library(dplyr)
library(grid)
library(gridExtra)
library(plotly)
library(xlsx)

setwd("~/OneDrive/University of Helsinki/HUMI/phage")

# Load in viral data after filtering iwth cutoff criteria 
phage_table <- read_csv("VP2_all_selection2.csv")

# separate by timepoint

phage_table_b4 <- phage_table[grep(".*B4", phage_table$sample_id), ]
phage_table_b5<- phage_table[grep(".*B5", phage_table$sample_id), ]
phage_table_b7<- phage_table[grep(".*B7", phage_table$sample_id), ]
phage_table_b9<- phage_table[grep(".*B9", phage_table$sample_id), ]
phage_table_M<- phage_table[grep(".*M", phage_table$sample_id), ]
phage_table_F<- phage_table[grep(".*F", phage_table$sample_id), ]

# remove huge table
rm(phage_table)

# Determine which technology returned information
phage_B4 <- phage_table_b4 %>% mutate(tools = ifelse((DVF & VS2),"Both",
                                               ifelse(DVF, "DVF", "VS2")))
phage_B5 <- phage_table_b5 %>% mutate(tools = ifelse((DVF & VS2),"Both",
                                               ifelse(DVF, "DVF", "VS2")))
phage_B7 <- phage_table_b7 %>% mutate(tools = ifelse((DVF & VS2),"Both",
                                               ifelse(DVF, "DVF", "VS2")))
phage_B9 <- phage_table_b9 %>% mutate(tools = ifelse((DVF & VS2),"Both",
                                               ifelse(DVF, "DVF", "VS2")))
phage_M <- phage_table_M %>% mutate(tools = ifelse((DVF & VS2),"Both",
                                             ifelse(DVF, "DVF", "VS2")))
phage_F <- phage_table_F %>% mutate(tools = ifelse((DVF & VS2),"Both",
                                             ifelse(DVF, "DVF", "VS2")))
# Set up data for plotting
tool_B4 <- phage_B4 %>% count(tools) %>%  mutate(time = c("B4","B4","B4"))
tool_B5 <- phage_B5 %>% count(tools) %>%  mutate(time = c("B5","B5","B5"))
tool_B7 <- phage_B7 %>% count(tools) %>%  mutate(time = c("B7","B7","B7"))
tool_B9 <- phage_B9 %>% count(tools) %>%  mutate(time = c("B9","B9","B9"))
tool_M <- phage_M %>% count(tools) %>%  mutate(time = c("M","M","M"))
tool_F <- phage_F %>% count(tools) %>%  mutate(time = c("F","F","F"))

# Plotting data
tool_graphing_data <- bind_rows(tool_B4, tool_B5, tool_B7, tool_B9,tool_M,tool_F)

# Plot
tool_graph <- ggplot(tool_graphing_data, aes(fill=tools, y=n, x=time)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Comparing Tool Results")
tool_graph

# CheckV quality 

checkV_B4 <- phage_B4 %>% count(checkv_quality) %>%  mutate(time = rep("B4",3))
checkV_B5 <- phage_B5 %>% count(checkv_quality) %>%  mutate(time = rep("B5",3))
checkV_B7 <- phage_B7 %>% count(checkv_quality) %>%  mutate(time = rep("B7",3))
checkV_B9 <- phage_B9 %>% count(checkv_quality) %>%  mutate(time = rep("B9",3))
checkV_M <- phage_M %>% count(checkv_quality) %>%  mutate(time = rep("M",4))
checkV_F <- phage_F %>% count(checkv_quality) %>%  mutate(time = rep("F",3))

checkV_graphing_data <- bind_rows(checkV_B4, checkV_B5, checkV_B7, checkV_B9, checkV_M, checkV_F)

checkV_graph <- ggplot(checkV_graphing_data, aes(fill=checkv_quality, y=n, x=time)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Checkv Quality Distribution")
checkV_graph

# Provirus 

provirus_B4 <- phage_B4 %>% count(provirus) %>%  mutate(time = rep("B4",2))
provirus_B5 <- phage_B5 %>% count(provirus) %>%  mutate(time = rep("B5",2))
provirus_B7 <- phage_B7 %>% count(provirus) %>%  mutate(time = rep("B7",2))
provirus_B9 <- phage_B9 %>% count(provirus) %>%  mutate(time = rep("B9",2))
provirus_M <- phage_M %>% count(provirus) %>%  mutate(time = rep("M",2))
provirus_F <- phage_F %>% count(provirus) %>%  mutate(time = rep("F",2))

provirus_graphing_data <- bind_rows(provirus_B4, provirus_B5, provirus_B7, provirus_B9, provirus_M, provirus_F)

provirus_graph <- ggplot(provirus_graphing_data, aes(fill=provirus, y=n, x=time)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Provirus Distribution")
provirus_graph

# Warning Table

warning_B4 <- phage_B4 %>% count(warnings, name = "B4")
warning_B5 <- phage_B5 %>% count(warnings, name = "B5")
warning_B7 <- phage_B7 %>% count(warnings, name = "B7") 
warning_B9 <- phage_B9 %>% count(warnings, name = "B9")
warning_M <- phage_M %>% count(warnings, name = "M")
warning_F <- phage_F %>% count(warnings, name = "F")

# Because vector warnings is of differing lengths for timepoint
max_length <- 8
warning_B7B7 <- warning_B7$B7
warning_FF <- warning_F$F
length(warning_B7B7) <- max_length                      
length(warning_FF) <- max_length

#cbind the three vectors together into a data frame
warning_data <- cbind(warning_B4$warnings, warning_B4$B4, warning_B5$B5, 
                      warning_B7B7, warning_B9$B9, warning_M$M, warning_FF)
colnames(warning_data) <- c("Warning","B4","B5","B7","B9","M","F")
warning_data[is.na(warning_data)] <- 0

# Contig Length 

getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Highlight peak in histogram
annotations_b4 <- data.frame(
  x = getmode(phage_B4$contig_length),
  y = 7500,
  label = "Peak: "
)

# Plot contig lengths
contig_phage_4 <- ggplot(phage_B4, aes(x=contig_length)) + 
  geom_histogram(color="black", fill="white", bins = 100) +
  xlim(0, 20000) +
  ggtitle("B4") +
  geom_text(data = annotations_b4, aes(x = x, y = y, label = paste(label, x)), size = 3)
contig_phage_4

annotations_b5 <- data.frame(
  x = getmode(phage_B5$contig_length),
  y = 15000,
  label = "Peak: "
)

contig_phage_5 <- ggplot(phage_B5, aes(x=contig_length)) + 
  geom_histogram(color="black", fill="white", bins = 100)  +
  xlim(0, 20000) +
  ggtitle("B5") +
  geom_text(data = annotations_b5, aes(x = x, y = y, label = paste(label, x)), size = 3)
contig_phage_5

annotations_b7 <- data.frame(
  x = getmode(phage_B7$contig_length),
  y = 27000,
  label = "Peak: "
)

contig_phage_7 <- ggplot(phage_B7, aes(x=contig_length)) + 
  geom_histogram(color="black", fill="white", bins = 100) +
  xlim(0, 20000) +
  ggtitle("B7") +
  geom_text(data = annotations_b7, aes(x = x, y = y, label = paste(label, x)), size = 3)
contig_phage_7

annotations_b9 <- data.frame(
  x = getmode(phage_B9$contig_length),
  y = 50000,
  label = "Peak: "
)

contig_phage_9 <- ggplot(phage_B9, aes(x=contig_length)) + 
  geom_histogram(color="black", fill="white", bins = 100) +
  xlim(0, 20000) +
  ggtitle("B9") +
  geom_text(data = annotations_b9, aes(x = x, y = y, label = paste(label, x)), size = 3) 
contig_phage_9

annotations_M <- data.frame(
  x = getmode(phage_M$contig_length),
  y = 57000,
  label = "Peak: "
)

contig_phage_M <- ggplot(phage_M, aes(x=contig_length)) + 
  geom_histogram(color="black", fill="white", bins = 100) +
  xlim(0, 20000) +
  ggtitle("M") +
  geom_text(data = annotations_M, aes(x = x, y = y, label = paste(label, x)), size = 3)
contig_phage_M

annotations_F <- data.frame(
  x = getmode(phage_F$contig_length),
  y = 30000,
  label = "Peak: "
)

contig_phage_F <- ggplot(phage_F, aes(x=contig_length)) + 
  geom_histogram(color="black", fill="white", bins = 100) +
  xlim(0, 20000) +
  ggtitle("F") +
  geom_text(data = annotations_F, aes(x = x, y = y, label = paste(label, x)), size = 3) 
contig_phage_F

## Splitting by sample

# Save files by time-point

write.table(phage_table_b4, file = "phage_B4_complete", sep = ",", row.names = FALSE)
write.table(phage_table_b5, file = "phage_B5_complete", sep = ",", row.names = FALSE)
write.table(phage_table_b7, file = "phage_B7_complete", sep = ",", row.names = FALSE)
write.table(phage_table_b9, file = "phage_B9_complete", sep = ",", row.names = FALSE)
write.table(phage_table_M, file = "phage_M_complete", sep = ",", row.names = FALSE)
write.table(phage_table_F, file = "phage_F_complete", sep = ",", row.names = FALSE)


