library(dplyr)
library(tidyverse)
library(devtools)
library(arsenal)

# Load in data
setwd("~/OneDrive/University of Helsinki/HUMI/phage/prodigal_prep")
viral_log <- read_csv("Catalogue_viralSeq_V1.03.23.csv")
viral_criteria <- read_csv("~/OneDrive/University of Helsinki/HUMI/phage/phage_criteria/VP2_all_selection2.csv")

# Group by sample
viral_log_sum <- viral_log %>% 
  group_by(sample_id) %>% 
  tally() 

# Group by sample
viral_criteria_sum <- viral_criteria %>% 
  group_by(sample_id) %>% 
  tally() 

# Join tables
viral_sums <- viral_log_sum %>% 
  left_join(viral_criteria_sum, by = "sample_id") %>% 
  rename("catalogue" = "n.x",
         "original" = "n.y")

# Filter out to keep only those samples that lost info
missing_data <- viral_sums %>% 
  filter(catalogue != original | xor(is.na(catalogue), is.na(original)))

# Detemine missing contigs 

contig_sample_group <- viral_criteria %>% 
  filter(grepl('Perhe1080_M|Perhe1101_F|Perhe21_B7|Perhe501_B5|
                Perhe507_B5|Perhe513_M|Perhe617_B7|Perhe69_B5', sample_id)) %>% 
  group_by(sample_id,contig_id) %>% 
  tally() 

viral_log_1080 <- viral_criteria %>% 
  filter(sample_id == "Perhe1080_M") %>% 
  filter(contig_id == "k141_145941")

viral_log_1101 <- viral_criteria %>% 
  filter(sample_id == "Perhe1101_F") %>% 
  filter(contig_id == "k141_170644")

viral_log_21 <- viral_criteria %>% 
  filter(sample_id == "Perhe21_B7") %>% 
  filter(contig_id == "k141_27098")

viral_log_501 <- viral_criteria %>% 
  filter(sample_id == "Perhe501_B5") %>% 
  filter(contig_id == "k141_15727")

viral_log_513 <- viral_criteria %>% 
  filter(sample_id == "Perhe513_M") %>% 
  filter(contig_id == "k141_303125")

viral_log_617 <- viral_criteria %>% 
  filter(sample_id == "Perhe617_B7") %>% 
  filter(contig_id == "k141_32951")

viral_log_617_2 <- viral_criteria %>% 
  filter(sample_id == "Perhe617_B7") %>% 
  filter(contig_id == "k141_43941")

viral_log_69 <- viral_criteria %>% 
  filter(sample_id == "Perhe69_B5") %>% 
  filter(contig_id == "k141_29733")

# Perhe507 

viral_criteria_507_contigs <- viral_criteria %>% 
  filter(sample_id == "Perhe507_B5") %>% 
  select(contig_id) 

viral_log_507_contigs <- viral_log %>% 
  filter(sample_id == "Perhe507_B5") %>% 
  select(contig_id)

# Found through duplicated(viral_criteria_507_contigs) the repeated contig is k141_104703

viral_criteria_507_contigs_2 <- viral_criteria %>% 
  filter(sample_id == "Perhe507_B5") %>% 
  filter(contig_id == "k141_104703")

  
