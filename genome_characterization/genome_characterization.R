library(dplyr)
library(tidyverse)
library(devtools)
library(ggtree)
library(bactaxR)

# Set working directory 

setwd("~/OneDrive/University of Helsinki/HUMI/Task2/genome_characterization")

# Load in files and organize 

fastANI <- read_tsv("FastANI_allvsAll.txt", col_names = c("File1","File2","ANI","ref1","ref2"))

map_database <- read_csv("Database_content_summary 9.41.20.csv")

fastANI <- fastANI[,1:3]
colnames(fastANI) <- c("Isolate","file_name","ANI")

# Keep only isolates of interest
 
fastANI <- fastANI[fastANI$ANI >= 99,]
# Remove the fact that highest ANI is one isolate to itself
fastANI <- fastANI[fastANI$Isolate != fastANI$file_name, ] 

# The Phocaeicola mapping names dont have the string "genomic" in it map database account for that here 

fastANI$file_name <- gsub("GCA_024758205.1_ASM2475820v1_genomic.fna", "GCA_024758205.1_ASM2475820v1.fna", fastANI$file_name)
fastANI$file_name <- gsub("GCA_020091465.1_ASM2009146v1_genomic.fna", "GCA_020091465.1_ASM2009146v1.fna", fastANI$file_name)

# Mapping 

map_databse_mapping <- select(map_database, `Organism Scientific Name` ,`Organism Qualifier`, file_name)

anno_fastANI <- fastANI %>% 
  left_join(
    map_databse_mapping,
    by = "file_name"
  )

# Dendrogram generation 

# Load in data
fastANI <- read_tsv("FastANI_allvsAll.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
fastANI_3 <- read_tsv("FastANI_allvsAll_Isolate3.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
fastANI_14 <- read_tsv("FastANI_allvsAll_Isolate14.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
fastANI_M1B2 <- read_tsv("FastANI_allvsAll_M1B2.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
fastANI_M1B3 <- read_tsv("FastANI_allvsAll_M1B3.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
fastANI_M1B5 <- read_tsv("FastANI_allvsAll_M1B5.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
fastANI_M2B2 <- read_tsv("FastANI_allvsAll_M2B2.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
fastANI_M2B3 <- read_tsv("FastANI_allvsAll_IsolateM2B3.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
fastANI_MAG <- read_tsv("FastANI_allvsAll_M1B3_MAG.txt", col_names = c("File1","File2","ANI","ref1","ref2"))

# Mapping 
map_database <- read_csv("Database_content_summary 9.41.20.csv")
Mapping_ref <- map_database %>% select(file_name, `Organism Scientific Name`,`Organism Qualifier`, `Assembly Name`) %>%
  rename("File2"="file_name")

Mapping_ref$`Organism Scientific Name` <- gsub("^(\\S*\\s+\\S+).*", "\\1", Mapping_ref$`Organism Scientific Name`)
Mapping_ref$tax <- paste(Mapping_ref$`Organism Scientific Name` , Mapping_ref$`Organism Qualifier`, Mapping_ref$`Assembly Name`)
Mapping_ref <- Mapping_ref[,c(1,5)]

# Tree annotation

# Mapped the two columns separately and had to adjust the naming of the column in the map_database to allow me to do so 
# The following code is to generate data that the dendrogram ANI function likes as an input to build dendrograms. Need n^2 number of rows 
# where n is the number of unique species 

colnames(Mapping_ref) <- c("File2","tax")

anno_fastANI <- left_join(fastANI, Mapping_ref, by = "File2")
anno_fastANI <- anno_fastANI %>% mutate(tax = ifelse(is.na(tax), File2, tax))

anno_fastANI <- anno_fastANI[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI <- left_join(anno_fastANI, Mapping_ref, by = "File1")
final_anno_fastANI <- final_anno_fastANI %>% mutate(tax.y = ifelse(is.na(tax.y), File1, tax.y))

final_anno_fastANI <- final_anno_fastANI[,c(4,2,3)]

final_anno_fastANI$tax.y <- gsub("IsolateM2B2eurofins_merged.fna", "Bacteroides thetaiotaomicron Isolate M2B2 ", final_anno_fastANI$tax.y)
final_anno_fastANI$tax.y <- gsub("Isolate14eurofins_merged.fna", "Bacteroides thetaiotaomicron Isolate 14 ", final_anno_fastANI$tax.y)
final_anno_fastANI$tax.x <- gsub("IsolateM2B2eurofins_merged.fna", "Bacteroides thetaiotaomicron Isolate M2B2 ", final_anno_fastANI$tax.x)
final_anno_fastANI$tax.x <- gsub("Isolate14eurofins_merged.fna", "Bacteroides thetaiotaomicron Isolate 14 ", final_anno_fastANI$tax.x)

final_anno_fastANI <- final_anno_fastANI[grep("thetaiotaomicron",final_anno_fastANI$tax.y),]
final_anno_fastANI <- final_anno_fastANI[grep("thetaiotaomicron",final_anno_fastANI$tax.x),]

write.table(final_anno_fastANI, "theta_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)

colnames(Mapping_ref) <- c("File2","tax")

# 3

anno_fastANI_3 <- left_join(fastANI_3, Mapping_ref, by = "File2")
anno_fastANI_3 <- anno_fastANI_3[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI_3 <- left_join(anno_fastANI_3, Mapping_ref, by = "File1")
final_anno_fastANI_3 <- final_anno_fastANI_3[,c(4,2,3)]
final_anno_fastANI_3$tax.x <- final_anno_fastANI_3$tax.x %>% replace_na('Bacteroides caccae Isolate 3')
final_anno_fastANI_3$tax.y <- final_anno_fastANI_3$tax.y %>% replace_na('Bacteroides caccae Isolate 3')

final_anno_fastANI_3 <- final_anno_fastANI_3[grep("caccae",final_anno_fastANI_3$tax.y),]
final_anno_fastANI_3 <- final_anno_fastANI_3[grep("caccae",final_anno_fastANI_3$tax.x),]

write.table(final_anno_fastANI_3, "3_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)
colnames(Mapping_ref) <- c("File2","tax")

# 14

anno_fastANI_14 <- left_join(fastANI_14, Mapping_ref, by = "File2")
anno_fastANI_14 <- anno_fastANI_14[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI_14 <- left_join(anno_fastANI_14, Mapping_ref, by = "File1")
final_anno_fastANI_14 <- final_anno_fastANI_14[,c(4,2,3)]
final_anno_fastANI_14$tax.x <- final_anno_fastANI_14$tax.x %>% replace_na('Bacteroides thetaiotaomicron Isolate 14')
final_anno_fastANI_14$tax.y <- final_anno_fastANI_14$tax.y %>% replace_na('Bacteroides thetaiotaomicron Isolate 14')

final_anno_fastANI_14 <- final_anno_fastANI_14[grep("thetaiotaomicron",final_anno_fastANI_14$tax.y),]
final_anno_fastANI_14 <- final_anno_fastANI_14[grep("thetaiotaomicron",final_anno_fastANI_14$tax.x),]

write.table(final_anno_fastANI_14, "14_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)
colnames(Mapping_ref) <- c("File2","tax")

# M1B2

anno_fastANI_M1B2 <- left_join(fastANI_M1B2, Mapping_ref, by = "File2")
anno_fastANI_M1B2 <- anno_fastANI_M1B2[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI_M1B2 <- left_join(anno_fastANI_M1B2, Mapping_ref, by = "File1")
final_anno_fastANI_M1B2 <- final_anno_fastANI_M1B2[,c(4,2,3)]
final_anno_fastANI_M1B2$tax.x <- final_anno_fastANI_M1B2$tax.x %>% replace_na('Bacteroides faecis Isolate M1B2')
final_anno_fastANI_M1B2$tax.y <- final_anno_fastANI_M1B2$tax.y %>% replace_na('Bacteroides faecis Isolate M1B2')

final_anno_fastANI_M1B2 <- final_anno_fastANI_M1B2[grep("faecis",final_anno_fastANI_M1B2$tax.y),]
final_anno_fastANI_M1B2 <- final_anno_fastANI_M1B2[grep("faecis",final_anno_fastANI_M1B2$tax.x),]

write.table(final_anno_fastANI_M1B2, "M1B2_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)
colnames(Mapping_ref) <- c("File2","tax")

# M1B3

anno_fastANI_M1B3 <- left_join(fastANI_M1B3, Mapping_ref, by = "File2")
anno_fastANI_M1B3 <- anno_fastANI_M1B3[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI_M1B3 <- left_join(anno_fastANI_M1B3, Mapping_ref, by = "File1")
final_anno_fastANI_M1B3 <- final_anno_fastANI_M1B3[,c(4,2,3)]
final_anno_fastANI_M1B3$tax.x <- final_anno_fastANI_M1B3$tax.x %>% replace_na('Bacteroides uniformis Isolate M1B3')
final_anno_fastANI_M1B3$tax.y <- final_anno_fastANI_M1B3$tax.y %>% replace_na('Bacteroides uniformis Isolate M1B3')

final_anno_fastANI_M1B3 <- final_anno_fastANI_M1B3[grep("uniformis",final_anno_fastANI_M1B3$tax.y),]
final_anno_fastANI_M1B3 <- final_anno_fastANI_M1B3[grep("uniformis",final_anno_fastANI_M1B3$tax.x),]

write.table(final_anno_fastANI_M1B3, "M1B3_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)
colnames(Mapping_ref) <- c("File2","tax")

# M1B5

anno_fastANI_M1B5 <- left_join(fastANI_M1B5, Mapping_ref, by = "File2")
anno_fastANI_M1B5 <- anno_fastANI_M1B5[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI_M1B5 <- left_join(anno_fastANI_M1B5, Mapping_ref, by = "File1")
final_anno_fastANI_M1B5 <- final_anno_fastANI_M1B5[,c(4,2,3)]
final_anno_fastANI_M1B5$tax.x <- final_anno_fastANI_M1B5$tax.x %>% replace_na('Bacteroides cellulosilyticus Isolate M1B5')
final_anno_fastANI_M1B5$tax.y <- final_anno_fastANI_M1B5$tax.y %>% replace_na('Bacteroides cellulosilyticus Isolate M1B5')

final_anno_fastANI_M1B5 <- final_anno_fastANI_M1B5[grep("cellulosilyticus",final_anno_fastANI_M1B5$tax.y),]
final_anno_fastANI_M1B5 <- final_anno_fastANI_M1B5[grep("cellulosilyticus",final_anno_fastANI_M1B5$tax.x),]

write.table(final_anno_fastANI_M1B5, "M1B5_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)
colnames(Mapping_ref) <- c("File2","tax")

# M2B2

anno_fastANI_M2B2 <- left_join(fastANI_M2B2, Mapping_ref, by = "File2")
anno_fastANI_M2B2 <- anno_fastANI_M2B2[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI_M2B2 <- left_join(anno_fastANI_M2B2, Mapping_ref, by = "File1")
final_anno_fastANI_M2B2 <- final_anno_fastANI_M2B2[,c(4,2,3)]
final_anno_fastANI_M2B2$tax.x <- final_anno_fastANI_M2B2$tax.x %>% replace_na('Bacteroides thetaiotaomicron Isolate M2B2')
final_anno_fastANI_M2B2$tax.y <- final_anno_fastANI_M2B2$tax.y %>% replace_na('Bacteroides thetaiotaomicron Isolate M2B2')

final_anno_fastANI_M2B2 <- final_anno_fastANI_M2B2[grep("thetaiotaomicron",final_anno_fastANI_M2B2$tax.y),]
final_anno_fastANI_M2B2 <- final_anno_fastANI_M2B2[grep("thetaiotaomicron",final_anno_fastANI_M2B2$tax.x),]

write.table(final_anno_fastANI_M2B2, "M2B2_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)
colnames(Mapping_ref) <- c("File2","tax")

# M2B3

fastANI_M2B3$File1 <- gsub("_genomic","",fastANI_M2B3$File1)
fastANI_M2B3$File2 <- gsub("_genomic","",fastANI_M2B3$File2)

# Mapping 

anno_fastANI_M2B3 <- left_join(fastANI_M2B3, Mapping_ref, by = "File2")
anno_fastANI_M2B3 <- anno_fastANI_M2B3[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI_M2B3 <- left_join(anno_fastANI_M2B3, Mapping_ref, by = "File1")
final_anno_fastANI_M2B3 <- final_anno_fastANI_M2B3[,c(4,2,3)]
final_anno_fastANI_M2B3$tax.x <- final_anno_fastANI_M2B3$tax.x %>% replace_na('Phocaeicola vulgatus Isolate M2B3')
final_anno_fastANI_M2B3$tax.y <- final_anno_fastANI_M2B3$tax.y %>% replace_na('Phocaeicola vulgatus Isolate M2B3')

write.table(final_anno_fastANI_M2B3, "M2B3_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)

# MAG 

colnames(Mapping_ref) <- c("File2","tax")

anno_fastANI_MAG <- left_join(fastANI_MAG, Mapping_ref, by = "File2")
anno_fastANI_MAG <- anno_fastANI_MAG %>% mutate(tax = ifelse(is.na(tax), File2, tax))

anno_fastANI_MAG <- anno_fastANI_MAG[,c(1,6,3)]

colnames(Mapping_ref) <- c("File1","tax")
final_anno_fastANI_MAG <- left_join(anno_fastANI_MAG, Mapping_ref, by = "File1")
final_anno_fastANI_MAG <- final_anno_fastANI_MAG %>% mutate(tax.y = ifelse(is.na(tax.y), File1, tax.y))

final_anno_fastANI_MAG <- final_anno_fastANI_MAG[,c(4,2,3)]

final_anno_fastANI_MAG$tax.y <- gsub("MAGs_B.Uniformis.fna", "Bacteroides uniformis MAG", final_anno_fastANI_MAG$tax.y)
final_anno_fastANI_MAG$tax.y <- gsub("IsolateM1B3eurofins_merged.fna", "Bacteroides uniformis Isolate M1B3", final_anno_fastANI_MAG$tax.y)
final_anno_fastANI_MAG$tax.x <- gsub("MAGs_B.Uniformis.fna", "Bacteroides uniformis MAG", final_anno_fastANI_MAG$tax.x)
final_anno_fastANI_MAG$tax.x <- gsub("IsolateM1B3eurofins_merged.fna", "Bacteroides uniformis Isolate M1B3", final_anno_fastANI_MAG$tax.x)

final_anno_fastANI_MAG <- final_anno_fastANI_MAG[grep("uniformis",final_anno_fastANI_MAG$tax.y),]
final_anno_fastANI_MAG <- final_anno_fastANI_MAG[grep("uniformis",final_anno_fastANI_MAG$tax.x),]

write.table(final_anno_fastANI_MAG, "unifirmis_MAG_ANI.txt", append = F, sep = "\t",
            row.names = F, col.names = F)

colnames(Mapping_ref) <- c("File2","tax")

# Create ANI object 
ani_3 <- read.ANI(file = "3_ANI.txt")
ani_14 <- read.ANI(file = "14_ANI.txt")
ani_M1B2 <- read.ANI(file = "M1B2_ANI.txt")
ani_M1B3 <- read.ANI(file = "M1B3_ANI.txt")
ani_M1B5 <- read.ANI(file = "M1B5_ANI.txt")
ani_M2B2 <- read.ANI(file = "M2B2_ANI.txt")
ani_M2B3 <- read.ANI(file = "M2B3_ANI.txt")
ani_MAG <- read.ANI(file = "unifirmis_MAG_ANI.txt")
ani_theta <- read.ANI(file = "theta_ANI.txt")

# Check data
summary(ani_3) 
summary(ani_14)
summary(ani_M1B2)
summary(ani_M1B3)
summary(ani_M1B5)
summary(ani_M2B2)
summary(ani_M2B3)
summary(ani_theta)
summary(ani_MAG)

# Create dendrogram 
dend_3 <- ANI.dendrogram(bactaxRObject = ani_3, ANI_threshold = 99, label_size = 0.25)
dend_14 <- ANI.dendrogram(bactaxRObject = ani_14, ANI_threshold = 99, label_size = 0.25)
dend_M1B2 <- ANI.dendrogram(bactaxRObject = ani_M1B2, ANI_threshold = 99, label_size = 0.25)
dend_M1B3 <- ANI.dendrogram(bactaxRObject = ani_M1B3, ANI_threshold = 99, label_size = 0.25)
dend_M1B5 <- ANI.dendrogram(bactaxRObject = ani_M1B5, ANI_threshold = 99, label_size = 0.25)
dend_M2B2 <- ANI.dendrogram(bactaxRObject = ani_M2B2, ANI_threshold = 99, label_size = 0.25)
dend_M2B3 <- ANI.dendrogram(bactaxRObject = ani_M2B3, ANI_threshold = 99, label_size = 0.25)
dend_theta <- ANI.dendrogram(bactaxRObject = ani_theta, ANI_threshold = 99, label_size = 0.25)
dend_MAG <- ANI.dendrogram(bactaxRObject = ani_MAG, ANI_threshold = 99, label_size = 0.25)


# Alise code

#fastANI <- read_tsv("FastANI_allvsAll.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
#fastANI_M2B3 <- read_tsv("FastANI_IsolateM2B3vsAll.txt", col_names = c("File1","File2","ANI","ref1","ref2"))
#fastANI_total <- rbind(fastANI, fastANI_M2B3)
#fastANI <- fastANI_total
#fastANI <- fastANI[fastANI$File1 != fastANI$File2, ]

#map_database <- read_csv("Database_content_summary.csv")

# Keep only isolates of interest
#fastANI_isolates <- fastANI %>% filter(str_detect(File1, "^Isolate")) %>%
  #filter(ANI>99)

# Mapping 
#Mapping_ref <- map_database %>% select(file_name, `Organism Scientific Name`) %>%
 #rename("File2"="file_name")

#anno_fastANI <- left_join(fastANI_isolates, Mapping_ref, by = "File2")

