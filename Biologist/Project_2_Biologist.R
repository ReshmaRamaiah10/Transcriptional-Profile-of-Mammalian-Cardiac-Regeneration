# Author: Allison Choy
# BF528 Project 2
# Role: Biologist

setwd('/projectnb2/bf528/users/tinman_2022/project_2_rest/')

# 7.1. Description: This section of code creates a line plot for FKPM values of representative 
# sarcomere, mitochondria, and cell cycle genes that are selected for their differential expression 
# as per O'Meara et al. in Fig. 1D.

#load libraries
library(tidyverse)
library(ggpubr)
library(patchwork)
library(pheatmap)

# read in FPKM tracking tables as a tibble
Ad_1<-tibble(read.table('/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking', header = TRUE))
Ad_2<-tibble(read.table('/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking', header = TRUE))
P0_1<-tibble(read.table('/projectnb/bf528/users/tinman_2022/project_2_rest/genes.fpkm_tracking_RK', header = TRUE))
P0_2<-tibble(read.table('/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking', header = TRUE))
P4_1<-tibble(read.table('/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking', header = TRUE))
P4_2<-tibble(read.table('/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking', header = TRUE))
P7_1<-tibble(read.table('/project/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking', header = TRUE))
P7_2<-tibble(read.table('/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking', header = TRUE))

# representative gene names from paper in its respective plot
#samples<-c(Ad_1, Ad_2, P0_2, P4_1, P4_2, P7_1, P7_2)
sarc_genes<-c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
mit_genes<-c('Mpc1', 'Prdx3', 'Acat1', 'Echs1', 'Slc25a11', 'Phyh')
cc_genes<-c('Cdc7', 'E2f8', 'Cdk7', 'Cdc26', 'Cdc6', 'E2f1', 'Cdc27', 'Bora', 'Cdc45', 'Rad51', 'Aurkb', 'Cdc23')

# Determines the mean of each timepoint, returns a joined tibble
mean_sample_FPKM<-function(sample_1, sample_2, gene_list){
  a<-sample_2 %>% 
    filter(gene_short_name %in% gene_list) %>%
    mutate(FPKM2 = FPKM) %>% 
    select(gene_short_name, FPKM2)
  sample_1 %>% 
    filter(gene_short_name %in% gene_list) %>% 
    left_join(a, "gene_short_name") %>% 
    select(gene_short_name, FPKM, FPKM2) %>% 
    group_by(gene_short_name) %>% 
    summarize(FPKM_mean = mean(FPKM:FPKM2)) %>% # Stdev = sd(FPKM:FPKM2)
    # took out sd due to issues downstream. Will also be omitted in plot: n of 2.
  return()
}

# Joins all FPKM for all timepoints into one tibble
sample_merge<-function(sample_1, sample_2, sample_3, sample_4){
  sample_1 %>% 
    #select(gene_short_name, FPKM_mean) %>% 
    left_join(sample_2, "gene_short_name") %>% 
    left_join(sample_3, "gene_short_name") %>% 
    left_join(sample_4, "gene_short_name") %>% 
  return()
}

# Sarcomere
# Creating tibble of FPKMs
Ad_src<-mean_sample_FPKM(Ad_1, Ad_2, sarc_genes)
P0_src<-mean_sample_FPKM(P0_1, P0_2, sarc_genes)
P4_src<-mean_sample_FPKM(P4_1, P4_2, sarc_genes)
P7_src<-mean_sample_FPKM(P7_1, P7_2, sarc_genes)

#combining all timepoints for each gene for plotting
sarc_FPKM<-sample_merge(P0_src, P4_src, P7_src, Ad_src)
colnames(sarc_FPKM)<-c('Gene', 'P0', 'P4', 'P7', 'Ad')
sarc_FPKM<-sarc_FPKM %>% 
  pivot_longer(cols = P0:Ad, names_to = 'timepoint', values_to = 'FKPM')
sarc_FPKM<-sarc_FPKM %>% 
  mutate(FPKM_factor = factor(sarc_FPKM$timepoint, levels = c('P0', 'P4', 'P7', 'Ad', ordered = TRUE)))

# Line Plot for Sarcomere
sarc_plot<-ggplot(sarc_FPKM, aes(x = fct_inorder(timepoint), 
                                 y = FKPM, group = Gene)) +
  geom_line(aes(color = Gene)) +
  geom_point(aes(shape = Gene)) +
  theme_classic() +
  scale_colour_manual(name='Gene',
                      labels = unique(sarc_FPKM$Gene),
                      values = c('red', 'orange', 'yellow', 'light green', 'green', 'cyan', 'magenta')) +
  scale_shape_manual(name = 'Gene',
                     labels = unique(sarc_FPKM$Gene),
                     values = c(1:7)) +
  labs(x = 'In vivo Maturation Timepoint', y = 'FPKM', title = 'Sarcomere') +
  theme(axis.title.x = element_text(vjust = -2), axis.text.x = element_text(angle = 45, vjust = -0.05),
        plot.title = element_text(hjust = 0.5), legend.position = 'right')
sarc_leg<-get_legend(sarc_plot)
#sarc_plot

# Mitchondria
# Tibble of FPKM
Ad_mit<-mean_sample_FPKM(Ad_1, Ad_2, mit_genes)
P0_mit<-mean_sample_FPKM(P0_2, P0_2, mit_genes)
P4_mit<-mean_sample_FPKM(P4_1, P4_2, mit_genes)
P7_mit<-mean_sample_FPKM(P7_1, P7_2, mit_genes)

# Mit tibble for plotting
mit_FPKM<-sample_merge(P0_mit, P4_mit, P7_mit, Ad_mit)
colnames(mit_FPKM)<-c('Gene', 'P0', 'P4', 'P7', 'Ad')
mit_FPKM<-mit_FPKM %>% 
  pivot_longer(cols = P0:Ad, names_to = 'timepoint', values_to = 'FKPM')
mit_FPKM<-mit_FPKM %>% 
  mutate(FPKM_factor = factor(mit_FPKM$timepoint, levels = c('P0', 'P4', 'P7', 'Ad', ordered = TRUE)))

# Line Plot for Mitochondria
mit_plot<-ggplot(mit_FPKM, aes(x = fct_inorder(timepoint), 
                                 y = FKPM, group = Gene)) +
  geom_line(aes(color = Gene)) +
  geom_point(aes(shape = Gene)) +
  theme_classic() +
  scale_colour_manual(name='Gene',
                      labels = unique(mit_FPKM$Gene),
                      values = c('red', 'orange', 'yellow', 'green', 'cyan', 'magenta')) +
  scale_shape_manual(name = 'Gene',
                     labels = unique(mit_FPKM$Gene),
                     values = c(1:5)) +
  labs(x = 'In vivo Maturation Timepoint', y = 'FPKM', title = 'Mitochondria') +
  theme(axis.title.x = element_text(vjust = -2), axis.text.x = element_text(angle = 45, vjust = -0.01),
        plot.title = element_text(hjust = 0.5), legend.position = 'right')
#mit_leg<-get_legend(mit_plot)
#mit_plot

# Cell Cycle
# CC Tibble
Ad_cc<-mean_sample_FPKM(Ad_1, Ad_2, cc_genes)
P0_cc<-mean_sample_FPKM(P0_2, P0_2, cc_genes)
P4_cc<-mean_sample_FPKM(P4_1, P4_2, cc_genes)
P7_cc<-mean_sample_FPKM(P7_1, P7_2, cc_genes)

# Cell Cycle tibble for plotting
cc_FPKM<-sample_merge(P0_cc, P4_cc, P7_cc, Ad_cc)
colnames(cc_FPKM)<-c('Gene', 'P0', 'P4', 'P7', 'Ad')
cc_FPKM<-cc_FPKM %>% 
  pivot_longer(cols = P0:Ad, names_to = 'timepoint', values_to = 'FKPM')
cc_FPKM<-cc_FPKM %>% 
  mutate(FPKM_factor = factor(cc_FPKM$timepoint, levels = c('P0', 'P4', 'P7', 'Ad', ordered = TRUE)))

# Line Plot for Cell Cycle
cc_plot<-ggplot(cc_FPKM, aes(x = fct_inorder(timepoint), 
                               y = FKPM, group = Gene)) +
  geom_line(aes(color = Gene)) +
  geom_point(aes(shape = Gene)) +
  theme_classic() +
  scale_colour_manual(name='Gene',
                      labels = unique(cc_FPKM$Gene),
                      values = c('red', 'orange', 'yellow', 'light green', 'green', 'cyan', 'sky blue', 'blue', 'purple', 'magenta', 'pink', 'black')) +
  scale_shape_manual(name = 'Gene',
                     labels = unique(cc_FPKM$Gene),
                     values = c(1:11)) +
  labs(x = 'In vivo Maturation Timepoint', y = 'FPKM', title = 'Cell Cycle') +
  theme(axis.title.x = element_text(vjust = -2), axis.text.x = element_text(angle = 45, vjust = -0.01),
        plot.title = element_text(hjust = 0.5))
#cc_plot

# Puts the three plots together - went with this option
ggarrange(sarc_plot, mit_plot, cc_plot, nrow = 2, ncol = 2, labels = 'AUTO')
ggsave('/projectnb/bf528/users/tinman_2022/project_2_rest/Fig1D_biologist.svg') #- saving it by hand has better dimensions
#sarc_plot|mit_plot|cc_plot - can be used but does not have label for parts A/B/C
#sarc_plot
#ggsave('/projectnb/bf528/users/tinman_2022/project_2_rest/sarc_FPKM.svg')
#mit_plot
#ggsave('/projectnb/bf528/users/tinman_2022/project_2_rest/mit_FPKM.svg')
#cc_plot
#ggsave('/projectnb/bf528/users/tinman_2022/project_2_rest/cc_FPKM.svg')
# 7.2. Description: Creates an augmented table from the DAVID analysis to include annotations that
# overlapped with those in the paper.

# read in tables as matrix

# initial attempt with tidy
#up_ref<-read_csv('paper_common_up.csv', col_names = TRUE, skip = 1)
#up<-read_delim('upregulated_clusters.txt', delim='\t', col_names = TRUE, skip = 1)
#down_ref<-read_csv('paper_common_down.csv', col_names = TRUE, skip = 1)
#down<-read_delim('downregulated_clusters.txt', delim='\t', col_names = TRUE, skip = 1)

up_ref2<-read.csv('paper_common_up.csv', header = FALSE)
names(up_ref2)<-up_ref2[2,]
up<-read.delim('upregulated_clusters.txt', header = FALSE, blank.lines.skip = TRUE)
names(up)<-up[2,]
down_ref<-read.csv('paper_common_down.csv', header = FALSE)
names(down_ref)<-down_ref[2,]
down<-read.delim('downregulated_clusters.txt', header = FALSE, blank.lines.skip = TRUE)
names(down)<-down[2,]

# convert them to tibble and clean up spacing/lines
up_ref<-up_ref2 %>% as_tibble() %>% 
  drop_na(Category) %>% 
  filter(Category != 'Category')

# tidy DAVID result and add Matching column to result
up<-up %>% as_tibble() %>% 
  drop_na(Category) %>% 
  filter(Category != 'Category') #%>% 
  #mutate(ifelse(up$Term %in% up_ref$Term,'Yes','No'))
#up_bool<-up$Term %in% up_ref$Term
up<-up %>% mutate('Match' = ifelse(up$Term %in% up_ref$Term, 'Yes', 'No'))
top_up_20<-up %>% filter(!str_detect(Category, 'Annotation')) %>% head(20) %>% select(2, 10:14)

# convert them to tibble and clean up spacing/lines
down_ref<-down_ref %>% as_tibble() %>% 
  drop_na(Category) %>% 
  filter(Category != 'Category')

# tidy DAVID result and add Matching column to result
down<-down %>% as_tibble() %>% 
  drop_na(Category) %>% 
  filter(Category != 'Category') #%>% 
  #mutate('Match' = ifelse(down$Term %in% down_ref$Term,'Yes','No'))
#down_bool<-down$Term %in% down_ref$Term
down<-down %>% mutate('Match'= ifelse(down$Term %in% down_ref$Term, 'Yes', 'No')) # works better this way
top_down_20<-down %>% filter(!str_detect(Category, 'Annotation')) %>% head(20) %>% select(2, 10:14)

# writes both into files
write.csv(up, file = 'upregulated ES.csv')
write.csv(top_up_20, file = 'Top 20 upregulated ES.csv')
write.csv(down, file = 'downregulated ES.csv')
write.csv(top_down_20, file = 'Top 20 downregulated ES.csv')

# 7.3. Description: Creates a FPKM matrix of all 8 samples concatenated from tracking files into a single dataframe.
# Then, subsetting the top 1000 genes that were differentially expressed from our P0 vs Ad analysis,
# creates a heatmap for comparison with Figure 2A.

# Tibble of samples with only tracking id, FPKM, and gene names
Ad_1_mat<-Ad_1 %>% 
  mutate(Ad_1_FPKM = FPKM) %>% 
  select(tracking_id, Ad_1_FPKM)
Ad_2_mat<-Ad_2 %>% 
  mutate(Ad_2_FPKM = FPKM) %>% 
  select(tracking_id, Ad_2_FPKM)
P0_1_mat<-P0_1 %>% 
  mutate(P0_1_FPKM = FPKM) %>% 
  select(tracking_id, gene_short_name, P0_1_FPKM)
P0_2_mat<-P0_2 %>% 
  mutate(P0_2_FPKM = FPKM) %>% 
  select(tracking_id, P0_2_FPKM)
P4_1_mat<-P4_1 %>% 
  mutate(P4_1_FPKM = FPKM) %>% 
  select(tracking_id, P4_1_FPKM)
P4_2_mat<-P4_2 %>% 
  mutate(P4_2_FPKM = FPKM) %>% 
  select(tracking_id, P4_2_FPKM)
P7_1_mat<-P7_1 %>% 
  mutate(P7_1_FPKM = FPKM) %>% 
  select(tracking_id, P7_1_FPKM)
P7_2_mat<-P7_2 %>% 
  mutate(P7_2_FPKM = FPKM) %>% 
  select(tracking_id, P7_2_FPKM)

# Joins all 8 into FPKM matrix, ranking them by tracking id
FPKM_matrix<-P0_1_mat %>% 
  distinct(tracking_id, .keep_all = TRUE) %>% 
  left_join(P0_2_mat, 'tracking_id') %>% 
  left_join(P4_1_mat, 'tracking_id') %>%
  left_join(P4_2_mat, 'tracking_id') %>% 
  left_join(P7_1_mat, 'tracking_id') %>% 
  left_join(P7_2_mat, 'tracking_id') %>% 
  left_join(Ad_1_mat, 'tracking_id') %>% 
  left_join(Ad_2_mat, 'tracking_id') %>% 
  arrange(rank(tracking_id)) # not sure if needed

# from the diff table of P0 vs Ad, obtains unique genes that are significant and ranked by q_value
cuffdiff<-tibble(read.table('/projectnb/bf528/users/tinman_2022/project_2_rest/gene_exp.diff', header = TRUE))
head(cuffdiff, 10)
cuffdiff_sorted<-cuffdiff %>% 
  filter(significant == 'yes') %>% 
  arrange(rank(q_value)) %>% 
  select(gene, log2.fold_change.:significant) %>% 
  distinct(gene, .keep_all = TRUE)

# subset of 1000 top genes that are differentially expressed
top_1000_cuffdiff<-cuffdiff_sorted[1:1000,]

# looks for any zeros in rows - not used
#FPKM_0<-apply(FPKM_matrix[3:ncol(FPKM_matrix),], 1, function(row) all(row != 0))
#FPKM_0<-FPKM_matrix %>% apply(MARGIN=1, function(row) all(row !=0 ))
#m<-FPKM_matrix[FPKM_0,]

# subsets FPKM matrix with the top 1000 genes and remove duplicates
FPKM_subset<-FPKM_matrix %>% 
  filter(rowSums(FPKM_matrix[,-c(1,2)])!=0) %>% 
  subset(gene_short_name %in% top_1000_cuffdiff$gene, select = c(2:10)) %>% 
  distinct(gene_short_name, .keep_all = TRUE)

# something I tried before
  #left_join(top_1000_cuffdiff, by=c('gene_short_name' = 'gene')) %>% 
  #filter(complete.cases(.)) %>% 
  #pivot_longer(cols = Ad_1_FPKM:P7_2_FPKM, names_to = 'Sample_timepoint', values_to = 'FPKM')
  
# subsetting only 200 genes out of 1000
FPKM_subset2<-FPKM_subset[1:200,]

# not sure if needed
genes<-pull(FPKM_subset[1])

# Plotting heatmap
p<-pheatmap(as.matrix(FPKM_subset2[-1]), color = colorRampPalette(c("blue", "brown", "yellow"))(100),
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
         labels_row = FPKM_subset2$gene_short_name, scale = 'row', fontsize_row = 1.8) #,
         #annotation_row = cluster_col)

# attempt at clustering
clust_mat <- as.matrix(FPKM_subset2)

library(dendextend)
#hclust for filtered_genes
clusters <- hclust(dist(FPKM_subset2[,2:9]), method = 'complete')
as.dendrogram(clusters) %>% plot(horiz = TRUE)
cluster_col<-cutree(tree = as.dendrogram(clusters), k = 13) %>% enframe(name = 'Sample', value = 'Cluster')
cluster_dendo<-tibble(cluster = case_when(cluster_col == 1 ~ 'cluster 1',
                                          cluster_col == 2 ~ 'cluster 2',
                                          cluster_col == 3 ~ 'cluster 3',
                                          cluster_col == 4 ~ 'cluster 4',
                                          cluster_col == 5 ~ 'cluster 5',
                                          cluster_col == 6 ~ 'cluster 6',
                                          cluster_col == 7 ~ 'cluster 7',
                                          cluster_col == 8 ~ 'cluster 8',
                                          cluster_col == 9 ~ 'cluster 9',
                                          cluster_col == 10 ~ 'cluster 10',
                                          cluster_col == 11 ~ 'cluster 11',
                                          cluster_col == 12 ~ 'cluster 12',
                                          cluster_col == 13 ~ 'cluster 13',
                                          TRUE ~ 'other'))

# Attempt at using geom_tile after pivot, but deleted some
# FPKM_plot<-FPKM_subset %>% 
#  ggplot(aes(x = Sample_timepoint, y = gene_short_name, fill = FPKM)) +
#  geom_tile() +
#  scale_fill_gradient(low="white", high="blue")
  #unique(FPKM_subset$gene_short_name)
#FPKM_plot

#FPKM_subset2<-FPKM_subset %>% ggplot() + geom_tile(aes(x = Sample_timepoint, y = gene_short_name, fill = FPKM))
#FPKM_subset2

# trying out other heatmaps
#pheatmap(as.matrix(FPKM_subset))
#heatmap(as.matrix(FPKM_subset), col = brewer.pal(15, RrBu))
