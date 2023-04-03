library(tidyverse)
library(dplyr)
library(scales)
library(plotly)
setwd("~/OneDrive/University of Helsinki/HUMI/Task2/AMR")

# Load in data and remove loose cutoff data
isolate3 <- read_tsv("isolate3.txt")
isolate3 <- isolate3[!grepl("Loose",isolate3$Cut_Off),]
isolate14 <- read_tsv("isolate14.txt")
isolate14 <- isolate14[!grepl("Loose",isolate14$Cut_Off),]
isolateM1B2 <- read_tsv("isolateM1B2.txt")
isolateM1B2 <- isolateM1B2[!grepl("Loose",isolateM1B2$Cut_Off),]
isolateM1B3 <- read_tsv("isolateM1B3.txt")
isolateM1B3 <- isolateM1B3[!grepl("Loose",isolateM1B3$Cut_Off),]
isolateM1B5 <- read_tsv("isolateM1B5.txt")
isolateM1B5 <- isolateM1B5[!grepl("Loose",isolateM1B5$Cut_Off),]
isolateM2B2 <- read_tsv("isolateM2B2.txt")
isolateM2B2 <- isolateM2B2[!grepl("Loose",isolateM2B2$Cut_Off),]
isolateM2B3 <- read_tsv("isolateM2B3.txt")
isolateM2B3 <- isolateM2B3[!grepl("Loose",isolateM2B3$Cut_Off),]

# List of the resistant drug classes for each isolate 
isolate3_cut <- data.frame(first_column  = c("glycopeptide antibiotic", "fluoroquinolone antibiotic","tetracycline antibiotic"))
colnames(isolate3_cut) <- c("Class")

isolate14_cut <- data.frame(first_column  = c("glycopeptide antibiotic", "fluoroquinolone antibiotic","tetracycline antibiotic"))
colnames(isolate14_cut) <- c("Class")

isolateM1B2_cut <- data.frame(first_column  = c("glycopeptide antibiotic", "fluoroquinolone antibiotic","tetracycline antibiotic"))
colnames(isolateM1B2_cut) <- c("Class")

isolateM1B3_cut <- data.frame(first_column  = c("glycopeptide antibiotic", "fluoroquinolone antibiotic","tetracycline antibiotic",
                                                 "cephalosporin","macrolide antibiotic"))
colnames(isolateM1B3_cut) <- c("Class")

isolateM1B5_cut <- data.frame(first_column  = c("glycopeptide antibiotic", "fluoroquinolone antibiotic","tetracycline antibiotic","macrolide antibiotic"))
colnames(isolateM1B5_cut) <- c("Class")

isolateM2B2_cut <- data.frame(first_column  = c("glycopeptide antibiotic", "fluoroquinolone antibiotic","tetracycline antibiotic"))
colnames(isolateM2B2_cut) <- c("Class")

isolateM2B3_cut <- data.frame (first_column  = c("glycopeptide antibiotic", "fluoroquinolone antibiotic","tetracycline antibiotic","macrolide antibiotic"))
colnames(isolateM2B3_cut) <- c("Class")

# Great dataframe to store data for plotting heatmap
heat_map <- data.frame(matrix(nrow = 35, ncol = 3))
colnames(heat_map) <- c("Isolate", "Dug_class", "TF")
heat_map$Isolate[1:5] <- "B. thetaiotaomicron (14)"
heat_map$Isolate[6:10] <- "B. cellulosilyticus (M1B5)"
heat_map$Isolate[11:15] <- "P. vulgatus (M2B3)"
heat_map$Isolate[16:20] <- "B. thetaiotaomicron (M2B2)"
heat_map$Isolate[21:25] <- "B. faecis (M1B2)"
heat_map$Isolate[26:30] <- "B. caccae (3)"
heat_map$Isolate[31:35] <- "B. uniformis (M1B3)"

# Add drug classes
drugs <- c("glycopeptide antibiotic", "fluoroquinolone antibiotic","tetracycline antibiotic",
           "cephalosporin","macrolide antibiotic")
heat_map$Dug_class[1:5] <- drugs
heat_map$Dug_class[6:10] <- drugs
heat_map$Dug_class[11:15] <- drugs
heat_map$Dug_class[16:20] <- drugs
heat_map$Dug_class[21:25] <- drugs
heat_map$Dug_class[26:30] <- drugs
heat_map$Dug_class[31:35] <- drugs

# Add drug class results
heat_map$TF[1:5] <- c(1,1,1,0,0)
heat_map$TF[6:10] <- c(1,1,1,0,0)
heat_map$TF[11:15] <- c(1,1,1,0,0)
heat_map$TF[16:20] <- c(1,1,1,1,1)
heat_map$TF[21:25] <- c(1,1,1,0,1)
heat_map$TF[26:30] <- c(1,1,1,0,0)
heat_map$TF[31:35] <- c(1,1,1,0,1)

# Plot
AMR_map <- ggplot(heat_map, aes(x=Dug_class, y=Isolate, fill= TF)) + 
  geom_tile() + scale_y_discrete(limits = c("B. uniformis (M1B3)", "B. caccae (3)","B. faecis (M1B2)","B. thetaiotaomicron (M2B2)","P. vulgatus (M2B3)","B. cellulosilyticus (M1B5)", "B. thetaiotaomicron (14)")) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) + ylab("Isolate") + xlab("Drug Class") + ggtitle("AMR Genes in Isolates") +
  scale_fill_gradient2(
    low = "black", 
    high = "red", 
    midpoint = .50)
AMR_map
