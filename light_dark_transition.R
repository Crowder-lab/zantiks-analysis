#!/usr/bin/env Rscript


# load library
source("utils.R")


# assay variable setup
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("light_dark_transition", suffixes)

# get files
main_files <- all_files[["main_files"]]
wildtype_files <- all_files[["wildtype_files"]]

# see if there's usable wildtype data
wildtype_exists <- length(wildtype_files) != 0
wildtype_only_names <- setdiff(names(wildtype_files), names(main_files))
wildtype_should_be_analyzed <- wildtype_exists && length(wildtype_only_names) > 0


# analysis function
analyze <- function(files) {
  analyzed_data <- list()
  for (prefix_name in names(files)) {
    group <- files[[prefix_name]]
    main_data <- load_zantiks(group$main, "dcidddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")
    genotypes <- load_genotypes(group$genotypes, group$fish_used, "across")

    long_data <- main_data %>%
      pivot_longer(
        cols = A1_Z1:A48_Z2,
        names_to = c("ARENA", "ZONE"),
        names_pattern = "A([0-9]+)_Z([1,2])",
        values_to = "DISTANCE"
      ) %>%
      mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
      relocate(ARENA, ZONE) %>%
      arrange(ARENA, ZONE, TIME)

    distance_data <- long_data %>%
      group_by(ARENA, BIN_NUM) %>%
      summarise(DISTANCE = sum(DISTANCE)) %>%
      ungroup() %>%
      mutate(minutes = as.integer(10 * BIN_NUM))

    analyzed_data[["distance"]][[prefix_name]] <- distance_data %>%
      attach_genotypes(genotypes)
  }

  analyzed_data
}


main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  prism_data <- xy_or_grouped_data(main_data[["distance"]][[prefix_name]], "DISTANCE", "minutes")
  write_csv(prism_data, file.path("data", "light_dark_transition", "output", paste0(prefix_name, "_DISTANCE.csv")))
}

# combine data (including wildtypes if possible)
if (wildtype_should_be_analyzed) {
  wildtype_data <- analyze(wildtype_files)
  all_names <- union(names(wildtype_files), names(main_files))

  combined_wildtype <- bind_rows(wildtype_data[["distance"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["distance"]], .id = "id")
  all_distance_data <- bind_rows(combined_wildtype, combined_main)
} else {
  all_names <- names(main_files)

  all_distance_data <- bind_rows(main_data[["distance"]], .id = "id")
}


# analyze and save combined data
for_prism <- add_numbering(all_distance_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "DISTANCE", "minutes")
write_csv(prism_data, file.path("data", "light_dark_transition", "output", "combined_DISTANCE.csv"))
