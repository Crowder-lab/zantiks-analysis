#!/usr/bin/env Rscript

# load library
source("utils.R")

# assay variable setup
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("sleep", suffixes)
main_files <- all_files[["main_files"]]
wildtype_files <- all_files[["wildtype_files"]]
wildtype_exists <- length(wildtype_files) != 0

# analysis function
analyze <- function(files) {
  analyzed_data <- list()
  for (prefix_name in names(files)) {
    # get group
    group <- files[[prefix_name]]

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
    analyzed_data[[prefix_name]] <- hourly_data %>%
      attach_genotypes(genotypes)
  }

  analyzed_data
}

# save main data files
main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  df <- main_data[[prefix_name]]
  prism_data <- df %>%
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
  write_csv(prism_data, file.path("data", "sleep", "output", paste0(prefix_name, ".csv")))
}

# handle wildtype data if it exists
if (wildtype_exists) {
  # we don't the wildtype data if it's also in our regular data
  wildtype_only_names <- setdiff(names(wildtype_files), names(main_files))
  wildtype_data <- analyze(wildtype_files)
  combined_wildtype <- bind_rows(wildtype_data, .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data, .id = "id")

  all_data <- bind_rows(combined_wildtype, combined_main)
  all_names <- union(names(wildtype_files), names(main_files))
} else {
  all_data <- bind_rows(main_data, .id = "id")
  all_names <- names(main_files)
}

# combine all data together
combined_data <- all_data %>%
  mutate(id_num = match(id, all_names)) %>% # convert each id (name from list) to int
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

# save the data
write_csv(combined_data, file.path("data", "sleep", "output", "combined.csv"))
