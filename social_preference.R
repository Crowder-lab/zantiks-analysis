#!/usr/bin/env Rscript


# load library
source("utils.R")


# assay variable setup
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("social_preference", suffixes)

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
    main_data <- load_zantiks(group$main, "didddddddddddddddddddddddddddddddddddddddddddddddddd")
    genotypes <- load_genotypes(group$genotypes, group$fish_used, "across")

    # make data longer form
    long_data <- main_data %>%
      pivot_longer(
        cols = T.A1.Z1:T.A10.Z5,
        names_to = c("ARENA", "ZONE"),
        names_pattern = "T.A([0-9]+).Z([1-5])",
        values_to = "SECONDS"
      ) %>%
      mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
      relocate(ARENA, ZONE) %>%
      arrange(ARENA, ZONE)

    social_preference_index_data <- long_data %>%
      group_by(ARENA) %>%
      summarise(
        social_preference_index = (
          sum(
            SECONDS[ZONE == 1] + # close to friend
              0.5 * SECONDS[ZONE == 2] -
              0.5 * SECONDS[ZONE == 4] -
              SECONDS[ZONE == 5] # far from friend
          ) /
            sum(SECONDS)
        )
      )

    zone_time_data <- long_data %>%
      group_by(ARENA, ZONE) %>%
      summarise(total_time = sum(SECONDS)) %>%
      ungroup()

    analyzed_data[["social preference index"]][[prefix_name]] <- social_preference_index_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["zone time"]][[prefix_name]] <- zone_time_data %>%
      attach_genotypes(genotypes)
  }

  analyzed_data
}


main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  # social preference index data
  prism_data <- column_data(main_data[["social preference index"]][[prefix_name]], "social_preference_index")
  write_csv(prism_data, file.path("data", "social_preference", "output", paste0(prefix_name, "_SOCIAL-PREFERENCE-INDEX.csv")))

  # zone time data
  prism_data <- xy_or_grouped_data(main_data[["zone time"]][[prefix_name]], "total_time", "ZONE")
  write_csv(prism_data, file.path("data", "social_preference", "output", paste0(prefix_name, "_ZONE-TIME.csv")))
}

# combine data (including wildtypes if possible)
if (wildtype_should_be_analyzed) {
  wildtype_data <- analyze(wildtype_files)
  all_names <- union(names(wildtype_files), names(main_files))

  # social preference index data
  combined_wildtype <- bind_rows(wildtype_data[["social preference index"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["social preference index"]], .id = "id")
  all_social_preference_index_data <- bind_rows(combined_wildtype, combined_main)

  # zone time data
  combined_wildtype <- bind_rows(wildtype_data[["zone time"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["zone time"]], .id = "id")
  all_zone_time_data <- bind_rows(combined_wildtype, combined_main)
} else {
  all_names <- names(main_files)

  # social preference index data
  all_social_preference_index_data <- bind_rows(main_data[["social preference index"]], .id = "id")

  # zone time data
  all_zone_time_data <- bind_rows(main_data[["zone time"]], .id = "id")
}


# analyze and save combined data
# social preference index
for_prism <- add_numbering(all_social_preference_index_data, all_names, "genotype")
prism_data <- column_data(for_prism, "social_preference_index")
write_csv(prism_data, file.path("data", "social_preference", "output", "combined_SOCIAL-PREFERENCE-INDEX.csv"))

# zone time
for_prism <- add_numbering(all_zone_time_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "total_time", "ZONE")
write_csv(prism_data, file.path("data", "social_preference", "output", "combined_ZONE-TIME.csv"))
