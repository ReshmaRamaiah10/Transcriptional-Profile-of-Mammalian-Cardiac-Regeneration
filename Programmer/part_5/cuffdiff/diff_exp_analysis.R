#load library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(xlsx)

#load in gene_diff data
gene_exp <- read.table("gene_exp.diff", header = TRUE)

#filter by significant(Y)
GEY <- filter(gene_exp, significant == "yes")


#write to xlsx to inspect GEY
write.xlsx(GEY, file = "significant_genes.xlsx",
           sheetName = "genes", append = FALSE)

#filter FPKM <1000 for both sample 1 and 2

#FILTER FPKM>1000
GEY_fpkm_sub_1000_1 <- filter(GEY, value_1 > 1000) 

GEY_fpkm_sub_1000_F <- filter(GEY_fpkm_sub_1000_1, value_2 > 1000) 

#write to xlsx to inspect GEY_fpkm_sub_1000_F
write.xlsx(GEY_fpkm_sub_1000_F, file = "DE_genes.xlsx",
           sheetName = "genes", append = FALSE)
