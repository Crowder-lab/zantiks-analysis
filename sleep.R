#!/usr/bin/env Rscript

# load library
source("utils.R")

# assay variable setup
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("sleep", suffixes)
main_files <- all_files[[1]]
wildtype_files <- all_files[[2]]

# analysis loop
all_data <- list()
for (i in seq_along(main_files)) {
  # get group
  group <- main_files[[i]]

  # get all requested data
  main_data <- load_zantiks(group$main, "dcidddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")
  genotypes <- load_genotypes(group$genotypes, group$fish_used, "down")

  long_data <- main_data %>%
    pivot_longer(
      cols = A1_Z1:A48_Z2,
      names_to = c("ARENA", "ZONE"),
      names_pattern = "A([0-9]+)_Z([1,2])",
      values_to = "DISTANCE"
    ) %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE)

  # bin numbers are better than seconds
  # seconds includes acclimation and tracking time :(
  hourly_data <- long_data %>%
    mutate(hour = as.integer(floor((BIN_NUM - 1) / 3300)) + 1) %>%
    group_by(ARENA, hour) %>%
    summarise(distance = sum(DISTANCE)) %>%
    ungroup()

  # add to list to combine later
  all_data[[names(main_files)[i]]] <- hourly_data %>%
    attach_genotypes(genotypes)

  prism_data <- hourly_data %>%
    attach_genotypes(genotypes) %>%
    arrange(ARENA, hour) %>%
    select(genotype, ARENA, hour, distance) %>%
    complete(genotype, ARENA = 1:200) %>%
    pivot_wider(
      names_from = c(genotype, ARENA),
      values_from = distance,
      names_glue = "{genotype}{ARENA}",
    ) %>%
    select(hour, paste0("WT", 1:200), paste0("HET", 1:200), paste0("HOM", 1:200)) %>%
    arrange(hour)

  # save the data!
  write_csv(prism_data, file.path("data", "sleep", "output", paste0(names(main_files)[i], ".csv")))
}

# combine all data together
combined_data <- bind_rows(all_data, .id = "id") %>%
  mutate(id_num = match(id, names(main_files))) %>% # convert each id (name from list) to int
  group_by(genotype) %>%
  arrange(id_num, ARENA) %>%
  mutate(unique_id = 96 * (id_num - 1) + ARENA - 1) %>% # both 1-indexed
  mutate(numbering = match(unique_id, sort(unique(unique_id)))) %>%
  ungroup() %>%
  select(-c(id, ARENA, id_num, unique_id)) %>%
  complete(genotype, numbering = 1:200) %>%
  pivot_wider(
    names_from = c(genotype, numbering),
    values_from = distance,
    names_glue = "{genotype}{numbering}"
  ) %>%
  select(hour, paste0("WT", 1:200), paste0("HET", 1:200), paste0("HOM", 1:200)) %>%
  arrange(hour)

write_csv(combined_data, file.path("data", "sleep", "output", "combined.csv"))
