#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
zantiks_indices <- as.vector(t(outer(LETTERS[1:6], sprintf("%02d", 1:8), FUN = paste0)))

df <- read_csv(args[1], show_col_types = FALSE) %>%
  pivot_longer(
    cols = -RUNTIME,
    names_to = c("coordinate", "ARENA"),
    names_pattern = "(X|Y)_A(\\d+)"
  ) %>%
  group_by(ARENA) %>%
  summarise(percent_missing_values = sum(is.na(value)) / n() * 100) %>%
  mutate(probably_wrong_zantiks = zantiks_indices[as.integer(ARENA)]) %>%
  select(c(ARENA, probably_wrong_zantiks, percent_missing_values)) %>%
  arrange(desc(percent_missing_values))

print(df, n = 48)
