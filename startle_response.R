#!/usr/bin/env Rscript

# load library
source("utils.R")

# set up variables for this assay
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt", xy = "_xy.csv")
all_files <- find_data("startle_response", suffixes)
main_files <- all_files[[1]]
wildtype_files <- all_files[[2]]

# analysis loop
for (i in seq_along(main_files)) {
  # get group
  group <- main_files[[i]]

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
      select(time_point, phase_id)

    # filter to ~2 s neighborhood around the times
    # make 30 fps examples match 5 fps examples
    # this is sort of stupid, only works for 30 and 5 fps
    # 5 fps ~= 50 rows, 30 fps ~= 300 rows
    # CHANGE IF YOU HAVE OTHER FRAMERATES
    if (nrow(xy_data) < 30000) { # 5 fps
      filtered_xy <- xy_data %>%
        cross_join(times) %>%
        filter(abs(RUNTIME - time_point) <= 1.4)
    } else { # 30 fps
      filtered_xy <- xy_data %>%
        filter(row_number() %% 6 == 0) %>% # decimate by a factor of 6
        cross_join(times) %>%
        filter(abs(RUNTIME - time_point) <= 1.4)
    }

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

    # add 200 ms bins relative to startle time
    # take only the first time value of each bin
    annotated_xy <- long_xy %>%
      mutate(relative_bin = as.integer(round((RUNTIME - time_point) * 5) * 200)) %>%
      group_by(ARENA, phase_id, relative_bin) %>%
      arrange(ARENA, phase_id, RUNTIME) %>%
      summarise(RUNTIME = first(RUNTIME), x = first(x), y = first(y)) %>%
      ungroup()

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
      arrange(ARENA, relative_bin)

    # prepare for graphpad prism
    prism_xy <- final_xy %>%
      select(genotype, ARENA, relative_bin, distance_step) %>%
      complete(genotype, ARENA = 1:200) %>%
      pivot_wider(
        names_from = c(genotype, ARENA),
        values_from = distance_step,
        names_glue = "{genotype}{ARENA}",
      ) %>%
      select(relative_bin, paste0("WT", 1:200), paste0("HET", 1:200), paste0("HOM", 1:200)) %>%
      arrange(relative_bin)

    # save the data!
    write_csv(prism_xy, file.path("data", "startle_response", "output", paste0(names(main_files)[i], "_", startle_type, ".csv")))
  }
}
