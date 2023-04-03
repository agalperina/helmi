library(dplyr)
library(tidyverse)
setwd("~/OneDrive/University of Helsinki/HUMI/Task2/prokka_annotations")

# Load prokka results

isolate_3 <- read_tsv("Isolate3_prokka.tsv")
isolate_14 <- read_tsv("Isolate14_prokka.tsv")
isolate_M1B2 <- read_tsv("IsolateM1B2_prokka.tsv")
isolate_M1B3 <- read_tsv("IsolateM1B3_prokka.tsv")
isolate_M1B5 <- read_tsv("IsolateM1B5_prokka.tsv")
isolate_M2B2 <- read_tsv("IsolateM2B2_prokka.tsv")
isolate_M2B3 <- read_tsv("IsolateM2B3_prokka.tsv")

# Determine number of annotations

hypo_genes_3 <- table(isolate_3$product)
# total annotations: 4519 
# hypothetical protein: 2687 
# percent hypothetical: 59.46
# percent annotated: 40.54

hypo_genes_14 <- table(isolate_14$product)
# total annotations: 4624
# hypothetical protein: 2605
# percent hypothetical: 56.34
# percent annotated: 43,66

hypo_genes_M1B2 <- table(isolate_M1B2$product)
# total annotations: 5394
# hypothetical protein: 3289
# percent hypothetical: 60.98
# percent annotated: 39.02

hypo_genes_M1B3 <- table(isolate_M1B3$product)
# total annotations: 4293
# hypothetical protein: 2442 
# percent hypothetical: 56.88
# percent annotated: 43.12

hypo_genes_M1B5 <- table(isolate_M1B5$product)
# total annotations: 5069
# hypothetical protein: 2856
# percent hypothetical: 56.34
# percent annotated: 43.66

hypo_genes_M2B2 <- table(isolate_M2B2$product)
# total annotations: 4771  
# hypothetical protein: 2709
# percent hypothetical: 56.78
# percent annotated: 43.22

hypo_genes_M2B3 <- table(isolate_M2B3$product)
# total annotations: 4150
# hypothetical protein: 2364
# percent hypothetical: 56.96
# percent annotated: 43.04
