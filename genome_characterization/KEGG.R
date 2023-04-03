library(tidyverse)
library(stringr)

setwd("~/OneDrive/University of Helsinki/HUMI/Task2/KEGG")

# Load data

isolate3 <- read_tsv("Isolate3_koala.txt")
isolate14 <- read_tsv("Isolate14_koala.txt")
isolateM1B2 <- read_tsv("IsolateM1B2_koala.txt")
isolateM1B3 <- read_tsv("IsolateM1B3_koala.txt")
isolateM1B5 <- read_tsv("IsolateM1B5_koala.txt")
isolateM2B2 <- read_tsv("IsolateM2B2_koala.txt")
isolateM2B3 <- read_tsv("IsolateM2B3_koala.txt")

KEGG_path <- read_csv("KEGG_pathwayko00540.csv", col_names = F)
colnames(KEGG_path) <- c("ID","Anno")

# Great dataframe to store data for plotting heatmap
heat_map <- data.frame(matrix(nrow = 434, ncol = 3)) 
colnames(heat_map) <- c("Isolate", "KEGG_ID", "TF")
heat_map$Isolate[1:62] <- "B. thetaiotaomicron (14)"
heat_map$Isolate[63:124] <- "B. cellulosilyticus (M1B5)"
heat_map$Isolate[125:186] <- "P. vulgatus (M2B3)"
heat_map$Isolate[187:248] <- "B. thetaiotaomicron (M2B2)"
heat_map$Isolate[249:310] <- "B. faecis (M1B2)"
heat_map$Isolate[311:372] <- "B. caccae (3)"
heat_map$Isolate[373:434] <- "B. uniformis (M1B3)"

# Add KEGG IDs
KEGG_ID <- KEGG_path$ID
heat_map$KEGG_ID[1:62] <- KEGG_ID
heat_map$KEGG_ID[63:124] <- KEGG_ID
heat_map$KEGG_ID[125:186] <- KEGG_ID
heat_map$KEGG_ID[187:248] <- KEGG_ID
heat_map$KEGG_ID[249:310] <- KEGG_ID
heat_map$KEGG_ID[311:372] <- KEGG_ID
heat_map$KEGG_ID[373:434] <- KEGG_ID

# Add KEGG results
heat_map$TF[1:62] <- ifelse(KEGG_path$ID %in% isolate14$KO, 1, 0)
heat_map$TF[63:124] <- ifelse(KEGG_path$ID %in% isolateM1B5$KO, 1, 0)
heat_map$TF[125:186] <- ifelse(KEGG_path$ID %in% isolateM2B3$KO, 1, 0)
heat_map$TF[187:248] <- ifelse(KEGG_path$ID %in% isolateM2B2$KO, 1, 0)
heat_map$TF[249:310] <- ifelse(KEGG_path$ID %in% isolateM1B2$KO, 1, 0)
heat_map$TF[311:372] <- ifelse(KEGG_path$ID %in% isolate3$KO, 1, 0)
heat_map$TF[373:434] <- ifelse(KEGG_path$ID %in% isolateM1B3$KO, 1, 0)

# Plot
ggplot(heat_map, aes(KEGG_ID, Isolate, fill= TF)) + 
  geom_tile() + scale_y_discrete(limits = c("B. uniformis (M1B3)", "B. caccae (3)","B. faecis (M1B2)","B. thetaiotaomicron (M2B2)","P. vulgatus (M2B3)","B. cellulosilyticus (M1B5)", "B. thetaiotaomicron (14)")) +
  theme(axis.text.x = element_text(angle = 90),plot.title = Yelement_text(hjust = 0.5)) + xlab("KEGG ID") + ylab("Isolate") + ggtitle("LPS Pathway Genes in Isolates") +
  scale_fill_gradient2(
    low = "black", 
    high = "red", 
    midpoint = .50)



