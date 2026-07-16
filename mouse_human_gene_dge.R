# Author: Kayla Mac
# Date: 7/13/26
# Data title: RNA sequencing analysis of hepatocellular carcinoma identified oxidative phosphorylation as a major pathologic feature
# Data set: GSE184733
# Purpose: find the mouse/human equivalents for each gene

pacman::p_load(tidyverse,
               ggplot2,
               ggrepel,
               edgeR,
               babelgene)

dge_human <- read_csv("Data/DGE/DGE.csv")

# get mouse equivalent of human genes in DGE file
orthologs <- orthologs(genes = dge_human$Gene, species = "mouse", human = TRUE)


dge_human_with_mouse <- dge_human |>
  left_join(orthologs |> # dge_human df join with orthologs, if no match in orthologs then NA
            select(human_symbol, symbol), # get human_symbol and symbol cols from orthologs
            by = c("Gene" = "human_symbol")) |> # match the Gene in dge_human to human_symbol in orthologs
  rename(mouse_gene = symbol) # rename the column to mouse_gene


