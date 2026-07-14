# Author: Kayla Mac
# Date: 7/13/26
# Data title: RNA sequencing analysis of hepatocellular carcinoma identified oxidative phosphorylation as a major pathologic feature
# Data set: GSE184733


#install.packages("pacman")

pacman::p_load(tidyverse,
               ggplot2,
               ggrepel,
               janitor,
               edgeR)


# function to compare whatever groups we want
compare_groups <- function(count_group, group, expt_name, ctrl_name) {
  
  group_factor <- factor(group, levels = c(ctrl_name, expt_name))
  
  y <- DGEList(counts = count_group, group = group_factor)
  
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  design <- model.matrix(~group_factor)
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, legacy = TRUE)
  
  qlf <- glmQLFTest(fit, coef = 2)
  
  results <- topTags(qlf, n = Inf)$table
  
  sorted_table <- results[order(results$logFC, decreasing = TRUE), ]
  
  return(sorted_table)
}


# read counts file
counts <- read_csv("Data/counts_cleaned.csv")


# compare nontumor and tumor
dge_groups <- factor(c(rep("Non-Tumor", 17), rep("Tumor", 17)))

dge <- compare_groups(counts, dge_groups, "Non-Tumor", "Tumor")

dge |>
  as.data.frame() |>
  write_csv("Data/DGE/DGE.csv")


