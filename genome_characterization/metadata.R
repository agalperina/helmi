library(plotly)
library(tidyverse)
require(gridExtra)

setwd("~/OneDrive/University of Helsinki/HUMI/Task2/metadata_info")
meta_phyloseq <- readRDS("PhySeq_BrakenHumGut_16.02.23.rds") 

# Relative abundance for each family to determine if sample from infant also exists in adult sample 

# Perhe 104

# Get sample
perhe_104 <- meta_phyloseq %>%
  subset_samples(meta_phyloseq@sam_data == "Perhe104") 
# Change object 
df.perhe_104 <- psmelt(perhe_104)
# Transform
df.perhe_104["Abundance"] <- transform(df.perhe_104["Abundance"], "compositional")
# Select only mother sample
df.perhe_104_M <- df.perhe_104[grepl("M",df.perhe_104$Sample_type),]

# Plot abundance 
comp_perhe_104_M <- ggplot(df.perhe_104_M, aes(x = Sample_type, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  ggtitle("M1 Genus") +
  scale_fill_discrete(breaks = "Bacteroides") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_104_M
#ggplotly(comp_perhe_104)

# Bacteroides abudance plot 
df.perhe_104_bact_M <- df.perhe_104_M[grepl("Bacteroides",df.perhe_104_M$genus),]
df.perhe_104_bact_M$isolate <- df.perhe_104_bact_M$species

# Highlight only species of interest 
for(i in 1:57){
  if(df.perhe_104_bact_M$species[i] == "Bacteroides faecis" | df.perhe_104_bact_M$species[i] == "Bacteroides uniformis" | df.perhe_104_bact_M$species[i] == "Bacteroides cellulosilyticus"){
    df.perhe_104_bact_M$isolate[i] = "interest"
  }
  else{
    df.perhe_104_bact_M$isolate[i] = "other"
  }
}

# Plot
legend_104 <- c("Bacteroides faecis", "Bacteroides uniformis", "Bacteroides cellulosilyticus")
comp_perhe_104_bact <- ggplot(df.perhe_104_bact_M, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  ggtitle("M1 Bacteroides") + 
  scale_color_manual(values = c(other = NA, interest = "black")) +
  scale_linetype_manual(values = c(other = "blank", interest = "dashed"), )+
  scale_fill_discrete(breaks = legend_104) + guides(color = "none", linetype = "none") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_104_bact

df.perhe_104_infant <- df.perhe_104[!grepl("M|F",df.perhe_104$Sample_type),]
df.perhe_104_infant$species[is.na(df.perhe_104_infant$species)] = "Unknown"
df.perhe_104_infant$isolate <- df.perhe_104_infant$species

for(i in 1:6260){
  if(df.perhe_104_infant$isolate[i] == "Bacteroides faecis" | df.perhe_104_infant$isolate[i] == "Bacteroides uniformis" | df.perhe_104_infant$isolate[i] == "Bacteroides cellulosilyticus"){
    next 
  }
  else{
    df.perhe_104_infant$isolate[i] = "other"
  }
}

# Species of interest in infants 
legend_104 <- c("Bacteroides faecis", "Bacteroides uniformis", "Bacteroides cellulosilyticus")
comp_perhe_104_bact_infant <- ggplot(df.perhe_104_infant, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  ggtitle("M1 Infant Bacteroides") + 
  scale_color_manual(values = c(other = NA, `Bacteroides faecis` = "black", `Bacteroides uniformis` = "red", `Bacteroides cellulosilyticus` = "blue")) +
  scale_linetype_manual(values = c(other = "blank", `Bacteroides faecis` = "solid", `Bacteroides uniformis` = "solid", `Bacteroides cellulosilyticus` = "solid"), )+
  scale_fill_discrete(breaks = legend_104) + guides(color = "none", linetype = "none") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_104_bact_infant

grid.arrange(comp_perhe_104, comp_perhe_104_bact, comp_perhe_104_bact_infant, ncol=3)

# Repeat for other samples
# Perhe 129

perhe_129 <- meta_phyloseq %>%
  subset_samples(meta_phyloseq@sam_data == "Perhe129") 
df.perhe_129 <- psmelt(perhe_129)
df.perhe_129["Abundance"] <- transform(df.perhe_129["Abundance"], "compositional")
df.perhe_129_M <- df.perhe_129[grepl("M",df.perhe_129$Sample_type),]

legend_129 <- c("Bacteroides", "Phocaeicola")
comp_perhe_129_M <- ggplot(df.perhe_129_M, aes(x = Sample_type, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  ggtitle("M2 Genus") +
  scale_fill_discrete(breaks = legend_129) + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_129
#ggplotly(comp_perhe_129)

df.perhe_129_bact_M <- df.perhe_129[grepl("Bacteroides",df.perhe_129$genus),]
df.perhe_129_bact_M$isolate <- df.perhe_129_bact_M$species

for(i in 1:57){
  if(df.perhe_129_bact_M$species[i] == "Bacteroides thetaiotaomicron"){
    df.perhe_129_bact_M$isolate[i] = "interest"
  }
  else{
    df.perhe_129_bact_M$isolate[i] = "other"
  }
}


comp_perhe_129_bact <- ggplot(df.perhe_129_bact_M, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  ggtitle("M2 Bacteroides") + 
  scale_color_manual(values = c(other = NA, interest = "black")) +
  scale_linetype_manual(values = c(other = "blank", interest = "dashed")) +
  scale_fill_discrete(breaks = "Bacteroides thetaiotaomicron") + guides(color = "none", linetype = "none") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_129_bact

df.perhe_129_pho_M <- df.perhe_129[grepl("Phocaeicola",df.perhe_129$genus),]
df.perhe_129_pho_M$isolate <- df.perhe_129_pho_M$species

for(i in 1:9){
  if(df.perhe_129_pho_M$species[i] == "Phocaeicola vulgatus"){
    df.perhe_129_pho_M$isolate[i] = "interest"
  }
  else{
    df.perhe_129_pho_M$isolate[i] = "other"
  }
}

comp_perhe_129_pho <- ggplot(df.perhe_129_pho_M, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  ggtitle("M2 Phocaeicola") + 
  scale_color_manual(values = c(other = NA, interest = "black")) +
  scale_linetype_manual(values = c(other = "blank", interest = "dashed")) +
  scale_fill_discrete(breaks = "Phocaeicola vulgatus") + guides(color = "none", linetype = "none") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_129_pho

df.perhe_129_infant <- df.perhe_129[!grepl("M|F",df.perhe_129$Sample_type),]
df.perhe_129_infant$species[is.na(df.perhe_129_infant$species)] = "Unknown"
df.perhe_129_infant$isolate <- df.perhe_129_infant$species

for(i in 1:6260){
  if(df.perhe_129_infant$isolate[i] == "Bacteroides thetaiotaomicron" | df.perhe_129_infant$isolate[i] == "Phocaeicola vulgatus"){
    next 
  }
  else{
    df.perhe_129_infant$isolate[i] = "other"
  }
}

legend_129 <- c("Bacteroides thetaiotaomicron", "Phocaeicola vulgatus")
comp_perhe_129_bact_infant <- ggplot(df.perhe_129_infant, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  ggtitle("M2 Longitude") + 
  scale_color_manual(values = c(other = NA, `Bacteroides thetaiotaomicron` = "black", `Phocaeicola vulgatus` = "red")) +
  scale_linetype_manual(values = c(other = "blank", `Bacteroides thetaiotaomicron` = "solid", `Phocaeicola vulgatus` = "solid"), )+
  scale_fill_discrete(breaks = legend_129) + guides(color = "none", linetype = "none") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_129_bact_infant

grid.arrange(comp_perhe_129, comp_perhe_129_bact, comp_perhe_129_pho,comp_perhe_129_bact_infant, ncol=2, nrow = 2)

# Perhe 25

perhe_25 <- meta_phyloseq %>%
  subset_samples(meta_phyloseq@sam_data == "Perhe25") 
df.perhe_25 <- psmelt(perhe_25)
df.perhe_25["Abundance"] <- transform(df.perhe_25["Abundance"], "compositional")
df.perhe_25_infant <- df.perhe_25[!grepl("M|F",df.perhe_25$Sample_type),]

comp_perhe_25_infant <- ggplot(df.perhe_25_infant, aes(x = Sample_type, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  ggtitle("14 Genus") +
  scale_fill_discrete(breaks = "Bacteroides") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_25_infant
#ggplotly(comp_perhe_25)

df.perhe_25_bact_infant <- df.perhe_25_infant[grepl("Bacteroides",df.perhe_25_infant$genus),]
df.perhe_25_bact_infant$isolate <- df.perhe_25_bact_infant$species

for(i in 1:285){
  if(df.perhe_25_bact_infant$species[i] == "Bacteroides thetaiotaomicron"){
    df.perhe_25_bact_infant$isolate[i] = "interest"
  }
  else{
    df.perhe_25_bact_infant$isolate[i] = "other"
  }
}

comp_perhe_25_bact_infant <- ggplot(df.perhe_25_bact_infant, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  ggtitle("14 Bacteroides") +
  scale_color_manual(values = c(other = NA, interest = "black")) +
  scale_linetype_manual(values = c(other = "blank", interest = "dashed")) +
  scale_fill_discrete(breaks = "Bacteroides thetaiotaomicron") + guides(color = "none", linetype = "none") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_25_bact_infant

df.perhe_25_m <- df.perhe_25[grepl("M",df.perhe_25$Sample_type),]
df.perhe_25_m <- df.perhe_25_m[grepl("Bacteroides",df.perhe_25_m$genus),]
df.perhe_25_m$isolate <- df.perhe_25_m$species

for(i in 1:57){
  if(df.perhe_25_m$species[i] == "Bacteroides thetaiotaomicron"){
    df.perhe_25_m$isolate[i] = "interest"
  }
  else{
    df.perhe_25_m$isolate[i] = "other"
  }
}

comp_perhe_25_bact_m <- ggplot(df.perhe_25_m, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  ggtitle("14 Longitude") + 
  scale_color_manual(values = c(other = NA, interst = "black")) +
  scale_linetype_manual(values = c(other = "blank", interest = "solid"), )+
  scale_fill_discrete(breaks = "Bacteroides thetaiotaomicron") + guides(color = "none", linetype = "none") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_25_bact_m

grid.arrange(comp_perhe_25_infant, comp_perhe_25_bact_infant,comp_perhe_25_bact_m, ncol=3)

# For Perhe 25 there is no B3 metadata even though the isolate came from timepoint B3

# Perhe 669

perhe_669 <- meta_phyloseq %>%
  subset_samples(meta_phyloseq@sam_data == "Perhe699") 
df.perhe_669 <- psmelt(perhe_669)
df.perhe_669["Abundance"] <- transform(df.perhe_669["Abundance"], "compositional")
df.perhe_669_infant<- df.perhe_669[!grepl("M|F",df.perhe_669$Sample_type),]

comp_perhe_669 <- ggplot(df.perhe_669_infant, aes(x = Sample_type, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  scale_fill_discrete(breaks = "Bacteroides") + theme(legend.position="bottom", legend.title = element_blank()) +
  ggtitle("3 Genus")
comp_perhe_669
ggplotly(comp_perhe_669)

df.perhe_669_bact_infant <- df.perhe_669_infant[grepl("Bacteroides",df.perhe_669_infant$genus),]
df.perhe_669_bact_infant$isolate <- df.perhe_669_bact_infant$species

for(i in 1:171){
  if(df.perhe_669_bact_infant$species[i] == "Bacteroides caccae"){
    df.perhe_669_bact_infant$isolate[i] = "interest"
  }
  else{
    df.perhe_669_bact_infant$isolate[i] = "other"
  }
}

comp_perhe_669_bact <- ggplot(df.perhe_669_bact_infant, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  scale_color_manual(values = c(other = NA, interest = "black")) +
  scale_linetype_manual(values = c(other = "blank", interest = "dashed")) +
  scale_fill_discrete(breaks = "Bacteroides caccae") + guides(color = "none", linetype = "none") +theme(legend.position="bottom", legend.title = element_blank()) +
  ggtitle("3 Bacteroides") 
comp_perhe_669_bact

df.perhe_669_m <- df.perhe_669[grepl("M",df.perhe_669$Sample_type),]
df.perhe_669_m <- df.perhe_669_m[grepl("Bacteroides",df.perhe_669_m$genus),]
df.perhe_669_m$species[is.na(df.perhe_669_m$species)] = "Unknown"
df.perhe_669_m$isolate <- df.perhe_669_m$species

for(i in 1:1565){
  if(df.perhe_669_m$isolate[i] == "Bacteroides caccae"){
    df.perhe_669_m$isolate[i] = "interest"
  } 
  else{
    df.perhe_669_m$isolate[i] = "other"
  }
}

comp_perhe_669_bact_m <- ggplot(df.perhe_669_m, aes(x = Sample_type, y = Abundance, fill = species, color = isolate, linetype = isolate)) + 
  geom_bar(stat = "identity") +
  ggtitle("3 Longitude") + 
  scale_color_manual(values = c(other = NA, interst = "black")) +
  scale_linetype_manual(values = c(other = "blank", interest = "solid"), )+
  scale_fill_discrete(breaks = "Bacteroides caccae") + guides(color = "none", linetype = "none") + theme(legend.position="bottom", legend.title = element_blank()) 
comp_perhe_669_bact_m

grid.arrange(comp_perhe_669, comp_perhe_669_bact, comp_perhe_669_bact_m, ncol=3)

