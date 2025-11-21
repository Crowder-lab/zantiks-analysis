#!/usr/bin/env Rscript
print("distance_traveled.R")
source("utils.R")


# assay variable setup
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("distance_traveled", suffixes)

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
    main_data <- load_zantiks(group$main, "ddicdddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")
    genotypes <- load_genotypes(group$genotypes, group$fish_used, "down")

    long_data <- main_data %>%
      pivot_longer(
        cols = A1_Z1:A48_Z2,
        names_to = c("ARENA", "ZONE"),
        names_pattern = "A([0-9]+)_Z([1,2])",
        values_to = "VALUE"
      ) %>%
      mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
      relocate(ARENA, ZONE) %>%
      arrange(ARENA, ZONE)

    distance_traveled_data <- long_data %>%
      group_by(ARENA) %>%
      summarise(distance_traveled = sum(VALUE[VARIABLE == "DISTANCE"]))

    percent_thigmotaxis_data <- long_data %>%
      group_by(ARENA) %>%
      filter(VARIABLE == "TIME") %>%
      summarise(percent_thigmotaxis = sum(VALUE[ZONE == 1]) / sum(VALUE) * 100)

    analyzed_data[["distance traveled"]][[prefix_name]] <- distance_traveled_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["percent thigmotaxis"]][[prefix_name]] <- percent_thigmotaxis_data %>%
      attach_genotypes(genotypes)
  }

  analyzed_data
}


main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  prism_data <- column_data(main_data[["distance traveled"]][[prefix_name]], "distance_traveled")
  write_csv(prism_data, file.path("data", "distance_traveled", "output", paste0(prefix_name, "_DISTANCE-TRAVELED.csv")))
  prism_data <- column_data(main_data[["percent thigmotaxis"]][[prefix_name]], "percent_thigmotaxis")
  write_csv(prism_data, file.path("data", "distance_traveled", "output", paste0(prefix_name, "_PERCENT-THIGMOTAXIS.csv")))
}

# combine data (including wildtypes if possible)
if (wildtype_should_be_analyzed) {
  wildtype_data <- analyze(wildtype_files)
  all_names <- union(names(wildtype_files), names(main_files))

  combined_wildtype <- bind_rows(wildtype_data[["distance traveled"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["distance traveled"]], .id = "id")
  all_distance_traveled_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["percent thigmotaxis"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["percent thigmotaxis"]], .id = "id")
  all_percent_thigmotaxis_data <- bind_rows(combined_wildtype, combined_main)
} else {
  all_names <- names(main_files)

  all_distance_traveled_data <- bind_rows(main_data[["distance traveled"]], .id = "id")
  all_percent_thigmotaxis_data <- bind_rows(main_data[["percent thigmotaxis"]], .id = "id")
}

# analyze and save combined data
for_prism <- add_numbering(all_distance_traveled_data, all_names, "genotype")
prism_data <- column_data(for_prism, "distance_traveled")
write_csv(prism_data, file.path("data", "distance_traveled", "output", "combined_DISTANCE-TRAVELED.csv"))

for_prism <- add_numbering(all_percent_thigmotaxis_data, all_names, "genotype")
prism_data <- column_data(for_prism, "percent_thigmotaxis")
write_csv(prism_data, file.path("data", "distance_traveled", "output", "combined_PERCENT-THIGMOTAXIS.csv"))
