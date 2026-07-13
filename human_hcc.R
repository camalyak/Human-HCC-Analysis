#install.packages("pacman")

pacman::p_load(tidyverse,
               ggplot2,
               ggrepel,
               janitor)

# read in counts file
counts <- read.delim("Data/GSE184733_counts_ALL.txt")

# rename X column to Gene
counts <- rename(counts, Gene = X)

