#!/usr/bin/env Rscript


# load library
source("utils.R")


# set up variables for this assay
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt", xy = "_xy.csv")
all_files <- find_data("startle_response", suffixes)

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
    main_data <- load_zantiks(group$main, "ddcdddddddddddddddddddddddddddddddddddddddddddddddd")
    xy_data <- load_xy(group$xy, 48)
    genotypes <- load_genotypes(group$genotypes, group$fish_used, "down")

    for (startle_type in c("STARTLE", "PREPULSE")) {
      # find out when the startles were
      times <- main_data %>%
        filter(PHASE == startle_type) %>%
        mutate(time_point = RUNTIME - 1) %>% # the RUNTIME gets recorded 1 s after the startle
        mutate(phase_id = paste0(startle_type, "_", 1:n())) %>%
        filter(phase_id %in% paste0(startle_type, "_", c(1, 2))) %>% # only first two startles and startle/prepulse
        select(time_point, phase_id)

      # filter to ~2 s neighborhood around the times
      filtered_xy <- xy_data %>%
        cross_join(times) %>%
        filter(abs(RUNTIME - time_point) <= 1.4)

      # make data long(er)-form: | RUNTIME | ARENA | time_point | phase_id | x | y |
      long_xy <- filtered_xy %>%
        pivot_longer(
          cols = -c(RUNTIME, time_point, phase_id),
          names_to = c("coordinate", "ARENA"),
          names_pattern = "(X|Y)_A(\\d+)"
        ) %>%
        mutate(ARENA = as.integer(ARENA)) %>%
        group_by(RUNTIME, ARENA, time_point, phase_id) %>%
        summarise(
          x = first(value),
          y = last(value),
        ) %>%
        ungroup()

      # add 100 ms bins relative to startle time
      # take only the first time value of each bin
      annotated_xy <- long_xy %>%
        mutate(relative_bin = as.integer(round((RUNTIME - time_point) * 10) * 100))# %>%

      # add distance in
      distance_xy <- annotated_xy %>%
        arrange(ARENA, phase_id, RUNTIME) %>%
        mutate(distance_step = sqrt((x - lag(x))^2 + (y - lag(y))^2))

      # collapse each phase into the mean
      averaged_xy <- distance_xy %>%
        group_by(ARENA, relative_bin) %>%
        summarise(distance_step = mean(distance_step, na.rm = TRUE)) %>%
        ungroup()

      # final filtering and column selection
      final_xy <- averaged_xy %>%
        filter(abs(relative_bin) <= 1000) %>%
        select(ARENA, relative_bin, distance_step) %>%
        attach_genotypes(genotypes) %>%
        group_by(ARENA) %>%
        # filter out fish that don't move at all around the startle from -500 ms to +500 ms
        filter(sum(distance_step[abs(relative_bin) <= 500]) != 0) %>%
        ungroup() %>%
        arrange(ARENA, relative_bin)

      # add to list to combine later
      analyzed_data[[startle_type]][[prefix_name]] <- final_xy
      analyzed_data[["ppi"]][[prefix_name]][[startle_type]] <- final_xy %>%
        group_by(ARENA) %>%
        summarise(genotype = first(genotype), !!sym(startle_type) := distance_step[relative_bin == 100])
    }

    startle_at_100ms <- analyzed_data[["ppi"]][[prefix_name]][["STARTLE"]]
    prepulse_at_100ms <- analyzed_data[["ppi"]][[prefix_name]][["PREPULSE"]]
    combined_at_100ms <- startle_at_100ms %>%
      inner_join(prepulse_at_100ms, by = join_by(ARENA, genotype)) %>%
      mutate(percent_ppi = (STARTLE - PREPULSE) / STARTLE * 100) %>%
      filter(percent_ppi <= 100 & percent_ppi >= 0) # also gets rid of infinites and NaNs
    analyzed_data[["ppi"]][[prefix_name]] <- combined_at_100ms
    }

  analyzed_data
}


# save main data files
main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  prism_data <- xy_or_grouped_data(main_data[["STARTLE"]][[prefix_name]], "distance_step", "relative_bin")
  write_csv(prism_data, file.path("data", "startle_response", "output", paste0(prefix_name, "_STARTLE.csv")))

  prism_data <- xy_or_grouped_data(main_data[["PREPULSE"]][[prefix_name]], "distance_step", "relative_bin")
  write_csv(prism_data, file.path("data", "startle_response", "output", paste0(prefix_name, "_PREPULSE.csv")))

  prism_data <- column_data(main_data[["ppi"]][[prefix_name]], "percent_ppi")
  write_csv(prism_data, file.path("data", "startle_response", "output", paste0(prefix_name, "_PERCENT-PPI.csv")))
}

# include wildtype as well
if (wildtype_should_be_analyzed) {
  wildtype_data <- analyze(wildtype_files)
  all_names <- union(names(wildtype_files), names(main_files))

  combined_wildtype <- bind_rows(wildtype_data[["STARTLE"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["STARTLE"]], .id = "id")
  all_startle_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["PREPULSE"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["PREPULSE"]], .id = "id")
  all_prepulse_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["ppi"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["ppi"]], .id = "id")
  all_percent_ppi_data <- bind_rows(combined_wildtype, combined_main)
} else {
  all_names <- names(main_files)

  all_startle_data <- bind_rows(main_data[["STARTLE"]], .id = "id")

  all_prepulse_data <- bind_rows(main_data[["PREPULSE"]], .id = "id")

  all_percent_ppi_data <- bind_rows(main_data[["ppi"]], .id = "id")
}

for_prism <- add_numbering(all_startle_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "distance_step", "relative_bin")
write_csv(prism_data, file.path("data", "startle_response", "output", "combined_STARTLE.csv"))

for_prism <- add_numbering(all_prepulse_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "distance_step", "relative_bin")
write_csv(prism_data, file.path("data", "startle_response", "output", "combined_PREPULSE.csv"))

for_prism <- add_numbering(all_percent_ppi_data, all_names, "genotype")
prism_data <- column_data(for_prism, "percent_ppi")
write_csv(prism_data, file.path("data", "startle_response", "output", "combined_PERCENT-PPI.csv"))
