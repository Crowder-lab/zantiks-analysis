#!/usr/bin/env Rscript


# load library
source("utils.R")


# set up variables for this assay
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("light_dark_preference_3wpf", suffixes)

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
    main_data <- load_zantiks(group$main, "dicdddddddddddddddddddddddd")
    genotypes <- load_genotypes(group$genotypes, group$fish_used, "across")

    # make data longer form
    long_data <- main_data %>%
      pivot_longer(
        cols = A1_Z1:A12_Z2,
        names_to = c("ARENA", "ZONE"),
        names_pattern = "A([0-9]+)_Z(1|2)",
        values_to = "VALUE"
      ) %>%
      mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
      relocate(ARENA, ZONE) %>%
      mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 1 ~ 3, .default = ZONE)) %>% # start swapping zone numbers
      mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 2 ~ 1, .default = ZONE)) %>% # now 1 is always dark
      mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 3 ~ 2, .default = ZONE)) %>% # and 2 is always light
      mutate(ZONE = c("dark", "light")[ZONE]) %>%
      attach_genotypes(genotypes) %>%
      arrange(ARENA, ZONE, TIME)

    # percent time spent in dark zone
    dark_time_data <- long_data %>%
      filter(ENDPOINT == "TIME_SPENT_IN_ZONE") %>%
      filter(ZONE == "dark") %>%
      mutate(percent_time = VALUE / 60 * 100) %>%
      mutate(minute = BIN_NUM) %>%
      select(genotype, ARENA, minute, percent_time)
    analyzed_data[["percent dark time"]][[prefix_name]] <- dark_time_data

    # light and dark distance
    distance_data <- long_data %>%
      filter(ENDPOINT == "DISTANCE_IN_ZONE") %>%
      group_by(ARENA, ZONE) %>%
      summarise(genotype = first(genotype), ZONE = first(ZONE), total_distance = sum(VALUE)) %>%
      ungroup()
    analyzed_data[["total distance"]][[prefix_name]] <- distance_data

    # light and dark time
    time_data <- long_data %>%
      filter(ENDPOINT == "TIME_SPENT_IN_ZONE") %>%
      group_by(ARENA, ZONE) %>%
      summarise(genotype = first(genotype), ZONE = first(ZONE), total_time = sum(VALUE)) %>%
      ungroup()
    analyzed_data[["total time"]][[prefix_name]] <- time_data
  }

  analyzed_data
}


# analyze and save each clutch of the main data
main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  # percent time spent in dark zone
  prism_data <- xy_or_grouped_data(main_data[["percent dark time"]][[prefix_name]], "percent_time", "minute")
  write_csv(prism_data, file.path("data", "light_dark_preference_3wpf", "output", paste0(prefix_name, "_PERCENT-DARK-TIME.csv")))

  # light and dark distance
  prism_data <- xy_or_grouped_data(main_data[["total distance"]][[prefix_name]], "total_distance", "ZONE")
  write_csv(prism_data, file.path("data", "light_dark_preference_3wpf", "output", paste0(prefix_name, "_TOTAL-DISTANCE.csv")))

  # light and dark time
  prism_data <- xy_or_grouped_data(main_data[["total time"]][[prefix_name]], "total_time", "ZONE")
  write_csv(prism_data, file.path("data", "light_dark_preference_3wpf", "output", paste0(prefix_name, "_TOTAL-TIME.csv")))
}


# combine data (including wildtypes if possible)
if (wildtype_exists) {
  wildtype_data <- analyze(wildtype_files)
  wildtype_only_names <- setdiff(names(wildtype_files), names(main_files))
  all_names <- union(names(wildtype_files), names(main_files))

  # percent time spent in dark zone
  combined_wildtype <- bind_rows(wildtype_data[["percent dark time"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["percent dark time"]], .id = "id")
  all_percent_dark_time_data <- bind_rows(combined_wildtype, combined_main)

  # light and dark distance
  combined_wildtype <- bind_rows(wildtype_data[["total distance"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["total distance"]], .id = "id")
  all_total_distance_data <- bind_rows(combined_wildtype, combined_main)

  # light and dark time
  combined_wildtype <- bind_rows(wildtype_data[["total time"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["total time"]], .id = "id")
  all_total_time_data <- bind_rows(combined_wildtype, combined_main)
} else {
  all_names <- names(main_files)

  # percent time spent in dark zone
  all_percent_dark_time_data <- bind_rows(main_data[["percent dark time"]], .id = "id")

  # light and dark distance
  all_total_distance_data <- bind_rows(main_data[["total distance"]], .id = "id")

  # light and dark time
  all_total_time_data <- bind_rows(main_data[["total time"]], .id = "id")
}


# analyze and save combined data
# percent time spent in dark zone
for_prism <- add_numbering(all_percent_dark_time_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "percent_time", "minute")
write_csv(prism_data, file.path("data", "light_dark_preference_3wpf", "output", "combined_PERCENT-DARK-TIME.csv"))

# light and dark distance
for_prism <- add_numbering(all_total_distance_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "total_distance", "ZONE")
write_csv(prism_data, file.path("data", "light_dark_preference_3wpf", "output", "combined_TOTAL-DISTANCE.csv"))

# light and dark time
for_prism <- add_numbering(all_total_time_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "total_time", "ZONE")
write_csv(prism_data, file.path("data", "light_dark_preference_3wpf", "output", "combined_TOTAL-TIME.csv"))
