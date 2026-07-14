# Author: Kayla Mac
# Date: 7/13/26
# Data title: RNA sequencing analysis of hepatocellular carcinoma identified oxidative phosphorylation as a major pathologic feature
# Data set: GSE184733


#install.packages("pacman")

pacman::p_load(tidyverse,
               ggplot2,
               ggrepel,
               janitor)


# read in tpm file
tpm <- read_csv("Data/tpm_matrix.csv")


# rename X column to Gene
tpm <- rename(tpm, Gene = ...1)


# convert to matrix
tpm_matrix <- tpm |>
  column_to_rownames(var = "Gene") |>
  as.matrix()


# scale the tpm matrix and normalize
tpm_matrix_scaled <- log2(tpm_matrix + 1)
tpm_matrix_scaled <- t(scale(t((tpm_matrix + 1))))


# put in ascending order and keep pairwise samples
tpm_matrix_scaled <- tpm_matrix_scaled[, order(
  as.numeric(sub("[NT]", "", colnames(tpm_matrix_scaled)))
)]


# convert back to df and save as csv
tpm_matrix_scaled |>
  as.data.frame() |> 
  rownames_to_column(var = "Gene") |>
  write_csv("Data/tpm_matrix_scaled.csv")


# read in counts file
counts <- read_delim("Data/GSE184733_counts_ALL.txt")


# rename X column to Gene
counts <- rename(counts, Gene = ...1)


# higher signal = MTARC genes, lower signal = MARCHF genes
fixed_counts <- counts |>
  mutate(row_mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) |>
  mutate(Gene = case_when(
    Gene == "1-Mar" & row_mean > 500 ~ "MTARC1",
    Gene == "2-Mar" & row_mean > 500 ~ "MTARC2",
    Gene == "1-Mar" & row_mean <= 500 ~ "MARCHF1",
    Gene == "2-Mar" & row_mean <= 500 ~ "MARCHF2",
    TRUE ~ Gene
  )) |>
  select(-row_mean)


# convert to matrix
counts_matrix <- fixed_counts |>
  column_to_rownames(var = "Gene") |>
  as.matrix()


# put in ascending order and keep pairwise samples
counts_matrix <- counts_matrix[, order(as.numeric(sub("[NT]", "", colnames(counts_matrix))))]


# save as csv file
counts_matrix |>
  as.data.frame() |>
  rownames_to_column(var = "Gene") |>
  write_csv("Data/counts_cleaned.csv")