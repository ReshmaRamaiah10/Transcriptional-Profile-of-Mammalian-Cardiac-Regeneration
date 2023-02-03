# Analyst
# Aneeq Husain

#importing libraries
library(tidyverse)

#opening the file
fp = '/projectnb2/bf528/users/tinman_2022/project_2_rest/gene_exp.diff'
gene_exp_df <- read.delim(fp)

#sort by q-value
gene_exp_df <- gene_exp_df[order(gene_exp_df$q_value),]

#subsetting the 10 most differentially expressed
top_ten_df <- gene_exp_df[1:10,c('gene','value_1','value_2','log2.fold_change.','p_value','q_value')]

#function to plot histograms
plot_hist <- function(df, title, color1){
  plot <- hist(df$log2.fold_change., breaks = 10, main = title, col = color1, xlab = 'Log2FC values', ylab = 'Count')
  return(plot)
}

#plot all genes
plot_hist(gene_exp_df, 'Histogram of Log2FC Values of all genes','coral2')

#df with significant genes
sig_df <- gene_exp_df %>% filter(significant == "yes")

# histogram of significant genes
plot_hist(sig_df, 'Histogram of Log2FC Values of significant genes','aquamarine4')

#separating up and down-regulated genes
upreg_df <- sig_df %>% filter(log2.fold_change. >0)
downref_df <- sig_df %>% filter(log2.fold_change. <0)

#write-out files with gene names
write(upreg_df$gene, "upregulated.txt")
write(downref_df$gene, "downregulated.txt")
