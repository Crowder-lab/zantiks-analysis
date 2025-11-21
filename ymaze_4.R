#!/usr/bin/env Rscript
print("ymaze_4.R")
# load library
source("utils.R")


# assay variables
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("ymaze_4", suffixes)

# get files
main_files <- all_files[["main_files"]]
wildtype_files <- all_files[["wildtype_files"]]

# see if there's usable wildtype data
wildtype_exists <- length(wildtype_files) != 0
wildtype_only_names <- setdiff(names(wildtype_files), names(main_files))
wildtype_should_be_analyzed <- wildtype_exists && length(wildtype_only_names) > 0


# analysis function
analyze <- function(files) {
  # set up a function to turn zone sequences into a turn direction
  zone_sequence_to_direction <- function(zone, next_zone) {
    case_when(
      zone == 1 & next_zone == 2 ~ "L",
      zone == 1 & next_zone == 3 ~ "R",
      zone == 2 & next_zone == 1 ~ "R",
      zone == 2 & next_zone == 3 ~ "L",
      zone == 3 & next_zone == 1 ~ "L",
      zone == 3 & next_zone == 2 ~ "R"
    )
  }

  analyzed_data <- list()
  for (prefix_name in names(files)) {
    group <- files[[prefix_name]]
    main_data <- load_zantiks(group$main, "dccici")
    genotypes <- load_genotypes(group$genotypes, group$fish_used, "across")

    zone_time_data <- main_data %>%
      arrange(ARENA, ZONE, TIME) %>%
      mutate(
        duration = ifelse(
          ACTION == "Enter_Zone" &
            lead(ACTION, 1) == "Exit_Zone" &
            !is.na(lead(ACTION, 1)) &
            lag(ACTION, 1) != "Enter_Zone",
          lead(TIME, 1) - TIME,
          0
        )
      ) %>%
      filter(duration > 0) %>%
      group_by(ARENA, ZONE) %>%
      summarise(zone_time = sum(duration, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(zone_type = ifelse(ZONE == 4, "Center", "Arms")) %>%
      group_by(ARENA, zone_type) %>%
      summarise(total_time = sum(zone_time)) %>%
      ungroup()

    turn_data <- main_data %>%
      filter(ZONE %in% c(1, 2, 3) & ACTION == "Enter_Zone") %>% # only entries into arms
      group_by(ARENA) %>%
      arrange(TIME) %>%
      mutate(next_zone = lead(ZONE), next_next_zone = lead(next_zone)) %>%
      filter(ZONE != next_zone) %>%
      mutate(turn_direction = zone_sequence_to_direction(ZONE, next_zone)) %>%
      mutate(next_turn = lead(turn_direction)) %>%
      mutate(next_next_turn = lead(next_turn)) %>%
      mutate(next_next_next_turn = lead(next_next_turn))

    spontaneous_alternation_percent_data <- turn_data %>%
      slice(-n()) %>%
      slice(-n()) %>%
      mutate(triad = paste0(ZONE, next_zone, next_next_zone)) %>%
      count(triad, name = "triad_count") %>%
      summarise(
        spontaneous_alternation_percent = sum(
          triad_count[triad %in% c("123", "231", "312", "321", "213", "132")]
        ) / sum(triad_count) * 100
      )

    tetragram_data <- turn_data %>%
      slice(-n()) %>%
      slice(-n()) %>%
      slice(-n()) %>%
      mutate(tetragram = paste0(
        turn_direction,
        next_turn,
        next_next_turn,
        next_next_next_turn
      )) %>%
      count(tetragram, name = "tetragram_count")

    turn_count_data <- turn_data %>%
      count(ZONE, name = "zone_count") %>%
      summarise(turn_count = sum(zone_count) - 1)

    alternation_percent_data <- tetragram_data %>%
      summarise(
        alternation_percent = sum(
          tetragram_count[tetragram %in% c("LRLR", "RLRL")]
        ) / sum(tetragram_count) * 100
      )

    repetition_percent_data <- tetragram_data %>%
      summarise(
        repetition_percent = sum(
          tetragram_count[tetragram %in% c("LLLL", "RRRR")]
        ) / sum(tetragram_count) * 100
      )

    analyzed_data[["zone time"]][[prefix_name]] <- zone_time_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["spontaneous alternation percent"]][[prefix_name]] <- spontaneous_alternation_percent_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["turn count"]][[prefix_name]] <- turn_count_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["alternation percent"]][[prefix_name]] <- alternation_percent_data %>%
      attach_genotypes(genotypes)
    analyzed_data[["repetition percent"]][[prefix_name]] <- repetition_percent_data %>%
      attach_genotypes(genotypes)
  }

  analyzed_data
}


main_data <- analyze(main_files)
for (prefix_name in names(main_files)) {
  prism_data <- xy_or_grouped_data(main_data[["zone time"]][[prefix_name]], "total_time", "zone_type")
  write_csv(prism_data, file.path("data", "ymaze_4", "output", paste0(prefix_name, "_ZONE-TIME.csv")))

  prism_data <- column_data(main_data[["spontaneous alternation percent"]][[prefix_name]], "spontaneous_alternation_percent")
  write_csv(prism_data, file.path("data", "ymaze_4", "output", paste0(prefix_name, "_SPONTANEOUS-ALTERNATION-PERCENT.csv")))

  prism_data <- column_data(main_data[["turn count"]][[prefix_name]], "turn_count")
  write_csv(prism_data, file.path("data", "ymaze_4", "output", paste0(prefix_name, "_TURN-COUNT.csv")))

  prism_data <- column_data(main_data[["alternation percent"]][[prefix_name]], "alternation_percent")
  write_csv(prism_data, file.path("data", "ymaze_4", "output", paste0(prefix_name, "_ALTERNATION-PERCENT.csv")))

  prism_data <- column_data(main_data[["repetition percent"]][[prefix_name]], "repetition_percent")
  write_csv(prism_data, file.path("data", "ymaze_4", "output", paste0(prefix_name, "_REPETITION-PERCENT.csv")))
}

if (wildtype_should_be_analyzed) {
  wildtype_data <- analyze(wildtype_files)
  all_names <- union(names(wildtype_files), names(main_files))

  combined_wildtype <- bind_rows(wildtype_data[["zone time"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["zone time"]], .id = "id")
  all_zone_time_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["spontaneous alternation percent"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["spontaneous alternation percent"]], .id = "id")
  all_spontaneous_alternation_percent_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["turn count"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["turn count"]], .id = "id")
  all_turn_count_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["alternation percent"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["alternation percent"]], .id = "id")
  all_alternation_percent_data <- bind_rows(combined_wildtype, combined_main)

  combined_wildtype <- bind_rows(wildtype_data[["repetition percent"]], .id = "id") %>%
    filter(genotype == "WT") %>%
    filter(id %in% wildtype_only_names)
  combined_main <- bind_rows(main_data[["repetition percent"]], .id = "id")
  all_repetition_percent_data <- bind_rows(combined_wildtype, combined_main)
} else {
  all_names <- names(main_files)

  all_zone_time_data <- bind_rows(main_data[["zone time"]], .id = "id")

  all_spontaneous_alternation_percent_data <- bind_rows(main_data[["spontaneous alternation percent"]], .id = "id")

  all_turn_count_data <- bind_rows(main_data[["turn count"]], .id = "id")

  all_alternation_percent_data <- bind_rows(main_data[["alternation percent"]], .id = "id")

  all_repetition_percent_data <- bind_rows(main_data[["repetition percent"]], .id = "id")
}

for_prism <- add_numbering(all_zone_time_data, all_names, "genotype")
prism_data <- xy_or_grouped_data(for_prism, "total_time", "zone_type")
write_csv(prism_data, file.path("data", "ymaze_4", "output", "combined_ZONE-TIME.csv"))

for_prism <- add_numbering(all_spontaneous_alternation_percent_data, all_names, "genotype")
prism_data <- column_data(for_prism, "spontaneous_alternation_percent")
write_csv(prism_data, file.path("data", "ymaze_4", "output", "combined_SPONTANEOUS-ALTERNATION-PERCENT.csv"))

for_prism <- add_numbering(all_turn_count_data, all_names, "genotype")
prism_data <- column_data(for_prism, "turn_count")
write_csv(prism_data, file.path("data", "ymaze_4", "output", "combined_TURN-COUNT.csv"))

for_prism <- add_numbering(all_alternation_percent_data, all_names, "genotype")
prism_data <- column_data(for_prism, "alternation_percent")
write_csv(prism_data, file.path("data", "ymaze_4", "output", "combined_ALTERNATION-PERCENT.csv"))

for_prism <- add_numbering(all_repetition_percent_data, all_names, "genotype")
prism_data <- column_data(for_prism, "repetition_percent")
write_csv(prism_data, file.path("data", "ymaze_4", "output", "combined_REPETITION-PERCENT.csv"))
