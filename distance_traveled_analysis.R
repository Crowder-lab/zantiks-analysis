#!/usr/bin/env Rscript

######## SETTING UP FUNCTIONS ########
library(tidyverse)

DEFAULT_GENOTYPES <- factor(c("WT", "HET", "HOM"), levels = c("WT", "HET", "HOM"))
PRISM_MAX_SEQUENCE <- 1:256


find_data <- function(assay_name, suffixes_needed) {
  # Validate inputs
  if (missing(assay_name) || is.null(assay_name) || assay_name == "") {
    stop("assay_name must be provided")
  }
  if (missing(suffixes_needed) || is.null(suffixes_needed)) {
    stop("suffixes_needed must be provided")
  }

  analyze_folder <- function(dir_name, suffixes_needed, is_main_folder) {
    unfiltered_csv_files <- list.files(dir_name, pattern = "\\.csv$", full.names = TRUE)
    # Validate that we have files
    if (is_main_folder && length(unfiltered_csv_files) == 0) {
      stop("No CSV files found in main directory")
    }

    # Remove other .csv files that are really from the same group
    csv_suffixes <- suffixes_needed[endsWith(unlist(suffixes_needed), ".csv")]
    for (csv_suffix in csv_suffixes) {
      unfiltered_csv_files <- unfiltered_csv_files[!endsWith(unfiltered_csv_files, csv_suffix)]
    }

    # Replace name after filtering
    main_files <- unfiltered_csv_files[!endsWith(unfiltered_csv_files, "_xy.csv")]
    # Extract prefixes from main data files
    prefixes <- tools::file_path_sans_ext(basename(main_files))

    # For each prefix, check if all suffix files exist
    file_groups <- list()
    for (prefix in prefixes) {
      current_group <- file_groups[[prefix]]
      current_group[["main"]] <- file.path(dir_name, paste0(prefix, ".csv"))

      for (suffix_name in names(suffixes_needed)) {
        suffix <- suffixes_needed[[suffix_name]]
        suffix_file <- file.path(dir_name, paste0(prefix, suffix))
        if (!file.exists(suffix_file)) {
          stop(paste0("Incomplete file group for prefix '", prefix, "': Missing file: ", suffix_file))
        }
        current_group[[suffix_name]] <- suffix_file
      }

      file_groups[[prefix]] <- current_group
    }

    file_groups
  }

  # Get data files from main directory
  main_dir <- file.path("raw zantiks data files", assay_name)
  if (!dir.exists(main_dir)) {
    stop(paste("Directory does not exist:", main_dir, ". Your current working directory is:", getwd()))
  }
  main_data <- analyze_folder(main_dir, suffixes_needed, TRUE)
  # Get wildtype data files
  wildtype_dir <- file.path("raw zantiks data files", assay_name, "wildtype")
  wildtype_data <- analyze_folder(wildtype_dir, suffixes_needed, FALSE)


  list(
    main_files = main_data,
    wildtype_files = wildtype_data
  )
}


load_zantiks <- function(file_path, col_types) {
  lines <- readLines(file_path)
  csv_text <- lines[4:(length(lines) - 1)]

  read_csv(I(csv_text), col_types = col_types)
}


load_xy <- function(file_path, num_arenas) {
  # X and Y for each arena plus the RUNTIME column
  read_csv(file_path, col_types = str_dup("d", num_arenas * 2 + 1))
}


load_genotypes <- function(genotyping_file, fish_used_file, counting_direction) {
  # read in data about which fish were used
  fish_used_data <- read.table(fish_used_file) %>%
    as.matrix()

  # create matrix to match well labels in genotype data
  label_matrix <- matrix(
    outer(LETTERS[1:8],
      sprintf("%02d", 1:12),
      FUN = paste0
    ),
    nrow = 8,
    ncol = 12
  )

  # get which wells were used
  wells_used <- label_matrix[fish_used_data == "x" | fish_used_data == "X"]

  # read in data
  genotype_data <- read_csv(genotyping_file) %>%
    # select(genotyping_well = Well, genotype = Cluster, clutch = Clutch) %>%
    select(genotyping_well = Well, genotype = Cluster) %>%
    mutate(row = str_extract(genotyping_well, "[A-H]"), column = as.integer(str_extract(genotyping_well, "[0-9]+"))) %>%
    mutate(genotyping_well = paste0(row, sprintf("%02d", column))) %>%
    filter(genotyping_well %in% wells_used)

  # sort by row or column
  if (counting_direction == "across") {
    genotype_data <- genotype_data %>%
      arrange(row, column)
  } else {
    genotype_data <- genotype_data %>%
      arrange(column, row)
  }

  # add row ids to join on
  genotype_data <- genotype_data %>%
    mutate(row_id = row_number())

  genotype_data
}


attach_genotypes <- function(data, genotypes) {
  attached_data <- data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    filter(!is.na(genotype)) %>% # unfilled row
    filter(genotype != "<Excluded>") %>% # HRM row that failed
    mutate(genotype = factor(genotype)) %>%
    select(genotype, names(data))

  # make genotypes be HOM, HET, WT if that seems right
  attached_genotypes <- unique(as.character(attached_data$genotype))
  if (all(attached_genotypes %in% unique(as.character(DEFAULT_GENOTYPES)))) {
    df <- attached_data %>%
      complete(genotype = DEFAULT_GENOTYPES) %>%
      mutate(genotype = fct_relevel(genotype, levels(DEFAULT_GENOTYPES)))
  } else {
    # otherwise leave them like they are
    df <- attached_data
  }

  df %>%
    arrange(genotype)
}


column_data <- function(data, values_column) {
  # get the genotypes that should be here
  genotype_levels <- levels(data$genotype)

  data_wide <- data %>%
    pivot_wider(
      id_cols = ARENA,
      names_from = genotype,
      values_from = !!sym(values_column),
    )

  # add genotype columnn if it doesn't exist
  for (col in genotype_levels) {
    if (!col %in% names(data_wide)) {
      data_wide[[col]] <- NA
    }
  }

  # get rid of invalid genotypes
  data_wide %>%
    select(all_of(genotype_levels))
}


xy_or_grouped_data <- function(data, values_column, non_genotype_column) {
  # get the genotypes present
  genotype_levels <- unique(as.character(levels(data$genotype)))

  data_wide <- data %>%
    # arenas now count from 1, no gaps
    group_by(genotype) %>%
    arrange(ARENA) %>%
    mutate(numbering = match(ARENA, sort(unique(ARENA)))) %>%
    ungroup() %>%
    select(-ARENA) %>%
    # make there be 256 arenas for padding
    complete(genotype, numbering = PRISM_MAX_SEQUENCE) %>%
    pivot_wider(
      names_from = c(genotype, numbering),
      values_from = !!sym(values_column),
      names_glue = "{genotype}_{numbering}",
    ) %>%
    select(
      !!sym(non_genotype_column),
      # select 256 columns of each genotype in order
      unlist(map(genotype_levels, ~ paste0(.x, "_", PRISM_MAX_SEQUENCE)))
    )

  data_wide
}


add_numbering <- function(data, id_values, grouping_column) {
  data %>%
    filter(!is.na(ARENA)) %>%
    mutate(id_num = match(id, id_values)) %>% # convert each id to int
    group_by(!!sym(grouping_column)) %>%
    arrange(id_num, ARENA) %>%
    mutate(unique_id = 96 * (id_num - 1) + ARENA - 1) %>% # both 1-indexed
    mutate(numbering = match(unique_id, sort(unique(unique_id)))) %>%
    ungroup() %>%
    select(-c(id, ARENA, id_num, unique_id)) %>%
    rename(ARENA = numbering)
}
######## SETTING UP FUNCTIONS ########


# assay variable setup
suffixes <- list(genotypes = "_genotypes.csv", fish_used = "_fish.txt")
all_files <- find_data("put your data here", suffixes)

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
  write_csv(prism_data, file.path("raw zantiks data files", "put your data here", "output", paste0(prefix_name, "_DISTANCE-TRAVELED.csv")))
  prism_data <- column_data(main_data[["percent thigmotaxis"]][[prefix_name]], "percent_thigmotaxis")
  write_csv(prism_data, file.path("raw zantiks data files", "put your data here", "output", paste0(prefix_name, "_PERCENT-THIGMOTAXIS.csv")))
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
write_csv(prism_data, file.path("raw zantiks data files", "put your data here", "output", "combined_DISTANCE-TRAVELED.csv"))

for_prism <- add_numbering(all_percent_thigmotaxis_data, all_names, "genotype")
prism_data <- column_data(for_prism, "percent_thigmotaxis")
write_csv(prism_data, file.path("raw zantiks data files", "put your data here", "output", "combined_PERCENT-THIGMOTAXIS.csv"))
