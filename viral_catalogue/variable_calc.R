library(dplyr)
library(tidyverse)
library(devtools)

setwd("~/OneDrive/University of Helsinki/HUMI/phage/technical_variables")

# Load in data
parent_vars <- read_csv("parental_technical_description_cohort_05.08.22.csv")
infant_vars <- read_csv("infant_technical_description_cohort_05.08.22.csv")
viral_data <- read_csv("~/OneDrive/University of Helsinki/HUMI/phage/phage_criteria/VP2_all_selection2.csv")

# Get viral counts
viral_counts <- viral_data %>% 
  group_by(sample_id) %>% 
  tally() 
colnames(viral_counts) <- c("Sample_ID", "Viral_Count")

# Sum of viral lengths
viral_length <- viral_data %>% 
  group_by(sample_id) %>% 
  tally(contig_length)
colnames(viral_length) <- c("Sample_ID", "Viral_Length")

# Split into parent and infant objects
viral_counts_parent <- viral_counts %>% 
  filter(grepl('M|F', Sample_ID)) %>% 
  arrange(match(Sample_ID, parent_vars$Sample_ID)) 
  
viral_counts_infants <- viral_counts %>% 
  filter(! grepl('M|F', Sample_ID)) %>% 
  arrange(match(Sample_ID, infant_vars$Sample_ID))

# 

viral_length_parent <- viral_length %>% 
  filter(grepl('M|F', Sample_ID)) %>% 
  arrange(match(Sample_ID, parent_vars$Sample_ID)) 

viral_length_infants <- viral_length %>% 
  filter(! grepl('M|F', Sample_ID)) %>% 
  arrange(match(Sample_ID, infant_vars$Sample_ID))

# Join tables
parent_vars <- parent_vars %>% 
  left_join(viral_counts_parent, by = "Sample_ID") %>% 
  left_join(viral_length_parent, by = "Sample_ID") %>% 
  drop_na(Viral_Count)

infant_vars <- infant_vars %>% 
  left_join(viral_counts_infants, by = "Sample_ID") %>% 
  left_join(viral_length_infants, by = "Sample_ID") %>% 
  drop_na(Viral_Count)

###  Numeric Variables -> Spearman ### 

## Parent ## 

# Home Storage

P_home_storage_count <- cor.test(parent_vars$Viral_Count, parent_vars$Storage_AtHome_Days, method = "spearman", exact = FALSE)
P_home_storage_length <- cor.test(parent_vars$Viral_Length, parent_vars$Storage_AtHome_Days, method = "spearman", exact = FALSE)

# Storage Lab

P_lab_storage_count <- cor.test(parent_vars$Viral_Count, parent_vars$Storage_Lab_Days, method = "spearman", exact = FALSE)
P_lab_storage_length <- cor.test(parent_vars$Viral_Length, parent_vars$Storage_Lab_Days, method = "spearman", exact = FALSE)

# Fecal Weight

P_fecal_weight_count <- cor.test(parent_vars$Viral_Count, parent_vars$extrac_fecalWeight, method = "spearman", exact = FALSE)
P_fecal_weight_length <- cor.test(parent_vars$Viral_Length, parent_vars$extrac_fecalWeight, method = "spearman", exact = FALSE)

# Extract Concentration

P_conc_count <- cor.test(parent_vars$Viral_Count, parent_vars$extrac_Conc, method = "spearman", exact = FALSE)
P_conc_length <- cor.test(parent_vars$Viral_Length, parent_vars$extrac_Conc, method = "spearman", exact = FALSE)

# Date

P_date_count <- cor.test(parent_vars$Viral_Count, parent_vars$Ext_Date_Unix, method = "spearman", exact = FALSE)
P_date_length <- cor.test(parent_vars$Viral_Length, parent_vars$Ext_Date_Unix, method = "spearman", exact = FALSE)

# Delivery 

P_delivery_count <- cor.test(parent_vars$Viral_Count, parent_vars$dist_delivery, method = "spearman", exact = FALSE)
P_delivery_length <- cor.test(parent_vars$Viral_Length, parent_vars$dist_delivery, method = "spearman", exact = FALSE)

# Combine 

parent_numeric <- list(P_conc, P_date, P_delivery, P_fecal_weight, P_lab_storage, P_home_storage)
names(parent_numeric) <- c("P_conc", "P_date","P_delivery", "P_fecal_weight", "P_lab_storage", "P_home_storage")

## Infant ## 

# Home Storage

I_home_storage_count <- cor.test(infant_vars$Viral_Count, infant_vars$Storage_AtHome_Days, method = "spearman", exact = FALSE)
I_home_storage_length <- cor.test(infant_vars$Viral_Length, infant_vars$Storage_AtHome_Days, method = "spearman", exact = FALSE)

# Storage Lab

I_lab_storage_count <- cor.test(infant_vars$Viral_Count, infant_vars$Storage_Lab_Days, method = "spearman", exact = FALSE)
I_lab_storage_length <- cor.test(infant_vars$Viral_Length, infant_vars$Storage_Lab_Days, method = "spearman", exact = FALSE)

# Fecal Weight

I_fecal_weight_count <- cor.test(infant_vars$Viral_Count, infant_vars$extrac_fecalWeight, method = "spearman", exact = FALSE)
I_fecal_weight_length <- cor.test(infant_vars$Viral_Length, infant_vars$extrac_fecalWeight, method = "spearman", exact = FALSE)

# Extract Concentration

I_conc_count <- cor.test(infant_vars$Viral_Count, infant_vars$extrac_Conc, method = "spearman", exact = FALSE)
I_conc_length <- cor.test(infant_vars$Viral_Length, infant_vars$extrac_Conc, method = "spearman", exact = FALSE)

# Date

I_date_count <- cor.test(infant_vars$Viral_Count, infant_vars$Ext_Date_Unix, method = "spearman", exact = FALSE)
I_date_length <- cor.test(infant_vars$Viral_Length, infant_vars$Ext_Date_Unix, method = "spearman", exact = FALSE)

# Combine 

infant_numeric <- list(I_conc, I_date, I_fecal_weight, I_lab_storage, I_home_storage)
names(infant_numeric) <- c("I_conc", "I_date", "I_fecal_weight", "I_lab_storage", "I_home_storage")

### Categorical variables -> Wilcoxin ###

## Parent ## 

# QC

P_QC_count <- pairwise.wilcox.test(parent_vars$Viral_Count, parent_vars$MGI_QC_Level, p.adjust.method = "BH")
P_QC_length <- pairwise.wilcox.test(parent_vars$Viral_Length, parent_vars$MGI_QC_Level, p.adjust.method = "BH")

# Plate

P_plate_count <- pairwise.wilcox.test(parent_vars$Viral_Count, parent_vars$extrac_PlateID, p.adjust.method = "BH")
P_plate_length <- pairwise.wilcox.test(parent_vars$Viral_Length, parent_vars$extrac_PlateID, p.adjust.method = "BH")

# Bristol

P_bristol_count <- pairwise.wilcox.test(parent_vars$Viral_Count, parent_vars$Ext_BristolScore, p.adjust.method = "BH")
P_bristol_length <- pairwise.wilcox.test(parent_vars$Viral_Length, parent_vars$Ext_BristolScore, p.adjust.method = "BH")

# Extractor

P_extractor_count <- pairwise.wilcox.test(parent_vars$Viral_Count, parent_vars$Ext_Extractor, p.adjust.method = "BH")
P_extractor_length <- pairwise.wilcox.test(parent_vars$Viral_Length, parent_vars$Ext_Extractor, p.adjust.method = "BH")

# Antibiotics

P_antibiotics_count <- pairwise.wilcox.test(parent_vars$Viral_Count, parent_vars$Coll_AntibioticTreatment, p.adjust.method = "BH")
P_antibiotics_length <- pairwise.wilcox.test(parent_vars$Viral_Length, parent_vars$Coll_AntibioticTreatment, p.adjust.method = "BH")

## Infant ## 

# QC

I_QC_count <- pairwise.wilcox.test(infant_vars$Viral_Count, infant_vars$MGI_QC_Level, p.adjust.method = "BH")
I_QC_length <- pairwise.wilcox.test(infant_vars$Viral_Length, infant_vars$MGI_QC_Level, p.adjust.method = "BH")

# Plate

I_Plate_count <- pairwise.wilcox.test(infant_vars$Viral_Count, infant_vars$extrac_PlateID, p.adjust.method = "BH")
I_Plate_length <- pairwise.wilcox.test(infant_vars$Viral_Length, infant_vars$extrac_PlateID, p.adjust.method = "BH")

# Bristol

I_Bristol_count <- pairwise.wilcox.test(infant_vars$Viral_Count, infant_vars$Ext_BristolScore, p.adjust.method = "BH")
I_Bristol_length <- pairwise.wilcox.test(infant_vars$Viral_Length, infant_vars$Ext_BristolScore, p.adjust.method = "BH")

# Extractor

I_Extractor_count <- pairwise.wilcox.test(infant_vars$Viral_Count, infant_vars$Ext_Extractor, p.adjust.method = "BH")
I_Extractor_length <- pairwise.wilcox.test(infant_vars$Viral_Length, infant_vars$Ext_Extractor, p.adjust.method = "BH")

