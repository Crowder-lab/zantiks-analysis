#!/usr/bin/env Rscript
print("sleep.R")
# load library
source("utils.R")


# assay variable setup
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("sleep", suffixes)

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

    distance_moved_data <- long_data %>%
      mutate(hour = as.integer(floor((BIN_NUM - 1) / 3300)) + 1) %>%
      filter(hour == 25) %>%
      group_by(ARENA) %>%
      summarise(distance_moved = sum(DISTANCE)) %>%
      ungroup()

    # this is slow to calculate
    thigmotaxis_data <- long_data %>%
      mutate(hour = as.integer(floor((BIN_NUM - 1) / 3300)) + 1) %>%
      filter(hour == 25) %>%
      group_by(ARENA, BIN_NUM) %>%
      summarise(primary_zone = ifelse(DISTANCE[ZONE == 1] > DISTANCE[ZONE == 2], 1, 2)) %>%
      summarise(percent_thigmotaxis = mean(primary_zone == 1) * 100)

    # add to list to combine later
    analyzed_data[["hourly"]][[prefix_name]] <- hourly_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["distance moved"]][[prefix_name]] <- distance_moved_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["percent thigmotaxis"]][[prefix_name]] <- thigmotaxis_data %>%
      attach_genotypes(genotypes)
  }

  analyzed_data
}


# analyze and save each clutch of the main data
main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  # hourly data
  prism_data <- xy_or_grouped_data(main_data[["hourly"]][[prefix_name]], "distance", "hour")
  write_csv(prism_data, file.path("data", "sleep", "output", paste0(prefix_name, "_HOURLY.csv")))

  # 1 hr distance moved
  prism_data <- column_data(main_data[["distance moved"]][[prefix_name]], "distance_moved")
  write_csv(prism_data, file.path("data", "sleep", "output", paste0(prefix_name, "_DISTANCE_MOVED.csv")))

  # percent thigmotaxis
  prism_data <- column_data(main_data[["percent thigmotaxis"]][[prefix_name]], "percent_thigmotaxis")
  write_csv(prism_data, file.path("data", "sleep", "output", paste0(prefix_name, "_PERCENT-THIGMOTAXIS.csv")))
}


# combine data (including wildtypes if possible)
if (wildtype_should_be_analyzed) {
  wildtype_data <- analyze(wildtype_files)
  all_names <- union(names(wildtype_files), names(main_files))

  combined_wildtype <- bind_rows(wildtype_data[["hourly"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["hourly"]], .id = "id")
  all_hourly_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["distance moved"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["distance moved"]], .id = "id")
  all_distance_moved_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["percent thigmotaxis"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["percent thigmotaxis"]], .id = "id")
  all_percent_thigmotaxis_data <- bind_rows(combined_wildtype, combined_main)
} else {
  all_names <- names(main_files)

  all_hourly_data <- bind_rows(main_data[["hourly"]], .id = "id")

  all_distance_moved_data <- bind_rows(main_data[["distance moved"]], .id = "id")

  all_percent_thigmotaxis_data <- bind_rows(main_data[["percent thigmotaxis"]], .id = "id")
}

for_prism <- add_numbering(all_hourly_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "distance", "hour")
write_csv(prism_data, file.path("data", "sleep", "output", "combined_HOURLY.csv"))

for_prism <- add_numbering(all_distance_moved_data, all_names, "genotype")
prism_data <- column_data(for_prism, "distance_moved")
write_csv(prism_data, file.path("data", "sleep", "output", "combined_DISTANCE-MOVED.csv"))

for_prism <- add_numbering(all_percent_thigmotaxis_data, all_names, "genotype")
prism_data <- column_data(for_prism, "percent_thigmotaxis")
write_csv(prism_data, file.path("data", "sleep", "output", "combined_PERCENT-THIGMOTAXIS.csv"))
