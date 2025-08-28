#!/usr/bin/env Rscript

# load library
source("utils.R")

# set up variables for this assay
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt", xy = "_xy.csv")
all_files <- find_data("startle_response", suffixes)
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
        group_by(ARENA) %>%
        # filter out fish that don't move at all around the startle from -400 ms to +400 ms
        filter(sum(distance_step[abs(relative_bin) <= 400]) != 0) %>%
        ungroup() %>%
        arrange(ARENA, relative_bin)

      # add to list to combine later
      analyzed_data[[startle_type]][[prefix_name]] <- final_xy
      analyzed_data[["ppi"]][[prefix_name]][[startle_type]] <- final_xy %>%
        group_by(ARENA) %>%
        summarise(genotype = first(genotype), !!sym(startle_type) := distance_step[relative_bin == 200])
    }

    startle_at_200ms <- analyzed_data[["ppi"]][[prefix_name]][["STARTLE"]]
    prepulse_at_200ms <- analyzed_data[["ppi"]][[prefix_name]][["PREPULSE"]]
    combined_at_200ms <- startle_at_200ms %>%
      inner_join(prepulse_at_200ms, by = join_by(ARENA, genotype)) %>%
      mutate(percent_ppi = (STARTLE - PREPULSE) / STARTLE * 100) %>%
      filter(percent_ppi <= 100 & percent_ppi >= 0) # also gets rid of infinites and NaNs
    analyzed_data[["ppi"]][[prefix_name]] <- combined_at_200ms
    }

  analyzed_data
}

# save main data files
main_data <- analyze(main_files)
for (startle_type in c("STARTLE", "PREPULSE")) {
  for (prefix_name in names(main_files)) {
    df <- main_data[[startle_type]][[prefix_name]]
    # prepare for graphpad prism
    prism_xy <- df %>%
      select(genotype, ARENA, relative_bin, distance_step) %>%
      complete(genotype, ARENA = 1:200) %>%
      pivot_wider(
        names_from = c(genotype, ARENA),
        values_from = distance_step,
        names_glue = "{genotype}{ARENA}",
      ) %>%
      select(relative_bin, paste0("WT", 1:200), paste0("HET", 1:200), paste0("HOM", 1:200)) %>%
      drop_na(relative_bin) %>%
      arrange(relative_bin)

    # save the data!
    write_csv(prism_xy, file.path("data", "startle_response", "output", paste0(prefix_name, "_", startle_type, ".csv")))
  }
}
for (prefix_name in names(main_files)) {
  df <- main_data[["ppi"]][[prefix_name]]
  prism_df <- df %>%
    select(ARENA, genotype, percent_ppi) %>%
    pivot_wider(names_from = genotype, values_from = percent_ppi) %>%
    select(WT, HET, HOM)

  # save the data
  write_csv(prism_df, file.path("data", "startle_response", "output", paste0(prefix_name, "_PERCENT-PPI.csv")))
}

# include wildtype as well
if (wildtype_exists) {
  wildtype_data <- analyze(wildtype_files)
}
for (startle_type in c("STARTLE", "PREPULSE")) {
  # if wildtype data exists
  if (wildtype_exists) {
    # remove main data from wildtype if present
    # also remove non-WT genotypes
    wildtype_only_names <- setdiff(names(wildtype_files), names(main_files))
    combined_wildtype <- bind_rows(wildtype_data[[startle_type]], .id = "id") %>%
      filter(genotype == "WT") %>%
      filter(id %in% wildtype_only_names)

    # get main data and combine
    combined_main <- bind_rows(main_data[[startle_type]], .id = "id")
    all_data <- bind_rows(combined_wildtype, combined_main)
    all_names <- union(names(wildtype_files), names(main_files))
  # if there's no wildtype data
  } else {
    all_data <- bind_rows(main_data[[startle_type]], .id = "id")
    all_names <- names(main_files)
  }

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
      values_from = distance_step,
      names_glue = "{genotype}{numbering}"
    ) %>%
    select(relative_bin, paste0("WT", 1:200), paste0("HET", 1:200), paste0("HOM", 1:200)) %>%
    drop_na(relative_bin) %>%
    arrange(relative_bin)

  write_csv(combined_data, file.path("data", "startle_response", "output", paste0("combined_", startle_type, ".csv")))
}
if (wildtype_exists) {
  combined_wildtype <- bind_rows(wildtype_data[["ppi"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["ppi"]], .id = "id")

  all_data <- bind_rows(combined_wildtype, combined_main)
  all_names <- union(names(wildtype_files), names(main_files))
} else {
  all_data <- bind_rows(main_data[["ppi"]], .id = "id")
  all_names <- names(main_files)
}
combined_data <- all_data %>%
  select(id, ARENA, genotype, percent_ppi) %>%
  pivot_wider(names_from = genotype, values_from = percent_ppi) %>%
  select(WT, HET, HOM)
write_csv(combined_data, file.path("data", "startle_response", "output", "combined_PERCENT-PPI.csv"))
