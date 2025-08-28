library(tidyverse)

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
    csv_suffixes = suffixes_needed[endsWith(unlist(suffixes_needed), ".csv")]
    for (csv_suffix in csv_suffixes) {
      unfiltered_csv_files <- unfiltered_csv_files[!endsWith(unfiltered_csv_files, csv_suffix)]
    }

    # Replace name after filtering
    main_files <- unfiltered_csv_files
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
  main_dir <- file.path("data", assay_name)
  if (!dir.exists(main_dir)) {
    stop(paste("Directory does not exist:", main_dir, ". Your current working directory is:", getwd()))
  }
  main_data <- analyze_folder(main_dir, suffixes_needed, TRUE)
  # Get wildtype data files
  wildtype_dir <- file.path("data", assay_name, "wildtype")
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
  data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    filter(genotype %in% c("WT", "HET", "HOM")) %>%
    arrange(genotype) %>%
    select(genotype, names(data))
}
