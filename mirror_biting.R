#!/usr/bin/env Rscript


# load library
source("utils.R")


# assay variable setup
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("mirror_biting", suffixes)

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
    main_data <- load_zantiks(group$main, "diddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiidddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")
    genotypes <- load_genotypes(group$genotypes, group$fish_used, "across")

    long_data <- main_data %>%
      pivot_longer(
        cols = D.A1.Z1:T.A20.Z3,
        names_to = c("CATEGORY", "ARENA", "ZONE"),
        names_pattern = "([D,C,T])\\.A([0-9]+)\\.Z([1-3])",
        values_to = "VALUE"
      ) %>%
      mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
      relocate(ARENA, ZONE) %>%
      arrange(ARENA, ZONE)

    mirror_distance_data <- long_data %>%
      filter(ZONE == 1 & CATEGORY == "D") %>%
      group_by(ARENA) %>%
      summarise(mirror_distance = sum(VALUE))

    mirror_time_data <- long_data %>%
      filter(ZONE == 1 & CATEGORY == "T") %>%
      group_by(ARENA) %>%
      summarise(mirror_time = sum(VALUE))

    analyzed_data[["mirror distance"]][[prefix_name]] <- mirror_distance_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["mirror time"]][[prefix_name]] <- mirror_time_data %>%
      attach_genotypes(genotypes)
  }

  analyzed_data
}


main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  prism_data <- column_data(main_data[["mirror distance"]][[prefix_name]], "mirror_distance")
  write_csv(prism_data, file.path("data", "mirror_biting", "output", paste0(prefix_name, "_MIRROR-DISTANCE.csv")))

  prism_data <- column_data(main_data[["mirror time"]][[prefix_name]], "mirror_time")
  write_csv(prism_data, file.path("data", "mirror_biting", "output", paste0(prefix_name, "_MIRROR-DISTANCE.csv")))
}

# combine data (including wildtypes if possible)
if (wildtype_should_be_analyzed) {
  wildtype_data <- analyze(wildtype_files)
  all_names <- union(names(wildtype_files), names(main_files))

  combined_wildtype <- bind_rows(wildtype_data[["mirror distance"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["mirror distance"]], .id = "id")
  all_mirror_distance_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["mirror time"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["mirror time"]], .id = "id")
  all_mirror_time_data <- bind_rows(combined_wildtype, combined_main)
} else {
  all_names <- names(main_files)

  all_mirror_distance_data <- bind_rows(main_data[["mirror distance"]], .id = "id")

  all_mirror_time_data <- bind_rows(main_data[["mirror time"]], .id = "id")
}


# analyze and save combined data
for_prism <- add_numbering(all_mirror_distance_data, all_names, "genotype")
prism_data <- column_data(for_prism, "mirror_distance")
write_csv(prism_data, file.path("data", "mirror_biting", "output", "combined_MIRROR-DISTANCE.csv"))

for_prism <- add_numbering(all_mirror_time_data, all_names, "genotype")
prism_data <- column_data(for_prism, "mirror_time")
write_csv(prism_data, file.path("data", "mirror_biting", "output", "combined_MIRROR-TIME.csv"))
