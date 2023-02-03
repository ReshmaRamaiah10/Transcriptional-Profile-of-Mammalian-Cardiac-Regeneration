#load library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(xlsx)

#load in fpkm data
genes.fpkm <- read.table("genes.fpkm_tracking_RK", header = TRUE)

#write to xlsx to inspect genes.fpkm
write.xlsx(genes.fpkm, file = "genes.fpkm_all.xlsx",
           sheetName = "genes.fpkm_all", append = FALSE)

#FILTER zeros
NZ <- filter(genes.fpkm, FPKM > 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001)

#write to xlsx to inspect
write.xlsx(NZ, file = "genes.fpkm_no_0.xlsx",
           sheetName = "genes.fpkm_no_0", append = FALSE)

#FPKM >1000 only
NZ_pro_1000 <-filter(NZ, FPKM >1000)

#histogram for all values
ggplot(genes.fpkm, aes(x = FPKM)) + geom_histogram(binwidth = 100000)

#hist No zero, pro 1000
ggplot(NZ_pro_1000, aes(x = FPKM)) + geom_histogram(binwidth = 35000)

#write to xlsx to inspect NZ_pro_1000
write.xlsx(NZ_pro_1000, file = "genes.fpkm_pro_1000.xlsx",
           sheetName = "genes.fpkm_pro_1000", append = FALSE)
