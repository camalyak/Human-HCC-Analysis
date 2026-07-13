#install.packages("pacman")

pacman::p_load(tidyverse,
               ggplot2,
               ggrepel,
               janitor)

# read in counts file
counts <- read.delim("Data/GSE184733_counts_ALL.txt")

counts <- rename(counts, Gene = X)
