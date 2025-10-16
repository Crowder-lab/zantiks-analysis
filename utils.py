import itertools
import os
import string

import cv2
import numpy as np
import pandas as pd
from numpy.typing import NDArray

DEFAULT_GENOTYPES = {"WT", "HET", "HOM"}
PRISM_MAX_SEQUENCE = range(1, 257)


def find_data(
    assay_name: str, suffixes_needed: dict[str, str]
) -> dict[str, dict[str, dict[str, str]]]:
    # TODO: add error conditions
    def analyze_folder(
        dir_name: str, suffixes_needed: dict[str, str], is_main_folder: bool
    ):
        unfiltered_csv_files = list(
            filter(lambda s: s.endswith(".csv"), os.listdir(dir_name))
        )
        # TODO: add error conditions
        csv_suffixes = tuple(
            filter(lambda s: s.endswith(".csv"), suffixes_needed.values())
        )
        unfiltered_csv_files = filter(
            lambda s: not s.endswith(csv_suffixes), unfiltered_csv_files
        )

        main_files = list(unfiltered_csv_files)
        prefixes = map(
            lambda s: os.path.splitext(s)[0], map(os.path.basename, main_files)
        )

        file_groups: dict[str, dict[str, str]] = {}
        for prefix in prefixes:
            file_groups[prefix] = {}
            file_groups[prefix]["main"] = os.path.join(dir_name, prefix + ".csv")

            for suffix_name in suffixes_needed.keys():
                suffix = suffixes_needed[suffix_name]
                suffix_file = os.path.join(dir_name, prefix + suffix)
                if not os.path.isfile(suffix_file):
                    # TODO: error
                    pass
                file_groups[prefix][suffix_name] = suffix_file

        return file_groups

    # Get data files from main directory
    main_dir = os.path.join("data", assay_name)
    if not os.path.isdir(main_dir):
        # TODO: error
        pass
    main_data = analyze_folder(main_dir, suffixes_needed, True)
    # Get wildtype data files
    wildtype_dir = os.path.join("data", assay_name, "wildtype")
    wildtype_data = analyze_folder(wildtype_dir, suffixes_needed, False)

    return {"main_files": main_data, "wildtype_files": wildtype_data}


def find_heatmap_data(assay_name: str, suffixes_needed: dict[str, str]):
    def analyze_folder(
        dir_name: str, suffixes_needed: dict[str, str], is_main_folder: bool
    ):
        main_files = list(filter(lambda s: s.endswith("xy.csv"), os.listdir(dir_name)))
        prefixes = map(
            lambda s: s.removesuffix("_xy.csv"), map(os.path.basename, main_files)
        )

        file_groups: dict[str, dict[str, str]] = {}
        for prefix in prefixes:
            file_groups[prefix] = {}
            file_groups[prefix]["main"] = os.path.join(dir_name, prefix + "_xy.csv")

            for suffix_name in suffixes_needed.keys():
                suffix = suffixes_needed[suffix_name]
                suffix_file = os.path.join(dir_name, prefix + suffix)
                if not os.path.isfile(suffix_file):
                    # TODO: error
                    pass
                file_groups[prefix][suffix_name] = suffix_file

        return file_groups

    # Get data files from main directory
    main_dir = os.path.join("data", assay_name)
    if not os.path.isdir(main_dir):
        # TODO: error
        pass
    main_data = analyze_folder(main_dir, suffixes_needed, True)
    # get wildtype data files
    wildtype_dir = os.path.join("data", assay_name, "wildtype")
    wildtype_data = analyze_folder(wildtype_dir, suffixes_needed, False)

    return {"main_files": main_data, "wildtype_files": wildtype_data}


def load_xy(file_path: str) -> pd.DataFrame:
    with open(file_path, "r") as f:
        df = pd.read_csv(f)
    return df


def load_genotypes(
    genotyping_file: str, fish_used_file: str, counting_direction: str
) -> pd.DataFrame:
    # read in data about which fish were used
    with open(fish_used_file, "r") as f:
        fish_used_data = np.asarray(
            list(map(lambda s: list(s.strip().replace(" ", "")), f.readlines()))
        )

    # create array to match well labels in genotype data
    label_array = np.asarray(
        list(
            map(
                lambda t: f"{t[0]}{t[1]:0>2}",
                itertools.product(string.ascii_uppercase[:8], range(1, 13)),
            )
        )
    ).reshape((8, 12))

    # get which wells were used
    wells_used = label_array[(fish_used_data == "x") | (fish_used_data == "X")]

    # read in data
    with open(genotyping_file, "r") as f:
        genotype_data = pd.read_csv(f)
    genotype_data = genotype_data[["Well", "Cluster"]].rename(
        columns={"Well": "genotyping_well", "Cluster": "genotype"}
    )
    genotype_data = genotype_data.assign(
        row=genotype_data["genotyping_well"].str.extract(r"([A-H])"),
        column=genotype_data["genotyping_well"].str.extract(r"([0-9]+)").astype(int),
    )
    genotype_data["genotyping_well"] = genotype_data["row"] + genotype_data[
        "column"
    ].astype(str).str.zfill(2)
    genotype_data = genotype_data.loc[genotype_data["genotyping_well"].isin(wells_used)]

    # sort by row or column
    if counting_direction == "across":
        genotype_data = genotype_data.sort_values(["row", "column"])
    else:
        genotype_data = genotype_data.sort_values(["column", "row"])

    # add row ids to join on
    genotype_data = genotype_data.reset_index(drop=True)

    return genotype_data


def load_arena_map(arena_map_file: str) -> NDArray:
    arena_map = cv2.imread(arena_map_file, cv2.IMREAD_COLOR_RGB)
    cleaned_arena_map = np.where(arena_map < 8, 0, arena_map)
    return cleaned_arena_map


def load_camera_params(camera_params_file: str) -> dict[str, NDArray]:
    camera_parameters = {}
    with open(camera_params_file, "rb") as f:
        temp_file = np.load(f)
        camera_parameters = {**temp_file}
    return camera_parameters


def get_arena_coords(arena_map: NDArray, arena: int) -> NDArray[bool]:
    # fmt: off
    color_cycle = np.asarray(((120, 120, 248),  # blue
                              (120, 248, 120),  # green
                              (120, 248, 248),  # cyan
                              (248, 120, 120),  # red
                              (248, 120, 248),  # pink
                              (248, 248, 120))) # yellow
    # fmt: on
    num_colors = 6
    darken_by = 8
    times_to_darken, color = divmod(arena - 1, num_colors)
    arena_color = color_cycle[color] - times_to_darken * darken_by
    return np.all(arena_map == arena_color, axis=-1)


def find_crop_coordinates(mask: NDArray[bool], buffer_scale: float = 0.1) -> tuple[int, int, int, int]:
    # find bounds of data
    row_min, row_max, col_min, col_max = None, None, None, None
    for i, row in enumerate(mask):
        row_sum = np.sum(row)
        if not row_min and row_sum > 0:
            row_min = i
        elif row_min and row_sum == 0:
            row_max = i
            break
    for i, col in enumerate(mask.T):
        col_sum = np.sum(col)
        if not col_min and col_sum > 0:
            col_min = i
        elif col_min and col_sum == 0:
            col_max = i
            break

    # add buffer to bounds for better viewing
    num_rows = row_max - row_min
    num_cols = col_max - col_min
    row_buffer_size = round(buffer_scale * num_rows)
    col_buffer_size = round(buffer_scale * num_cols)
    row_min -= row_buffer_size
    row_max += row_buffer_size
    col_min -= col_buffer_size
    col_max += col_buffer_size

    # error if we go out of bounds
    if row_min < 0 or col_min < 0 or row_max > mask.shape[0] or col_max > mask.shape[1]:
        raise ValueError("Buffer scale too big. Image bounds exceeded.")

    return row_min, row_max, col_min, col_max


def attach_genotypes(data: pd.DataFrame, genotypes: pd.DataFrame) -> pd.DataFrame:
    attached_data = data.join(genotypes, on="ARENA", how="left")
    attached_data = attached_data.dropna(subset="genotype")
    attached_data = attached_data.loc[attached_data["genotype"] != "<Excluded>"]
    g = ["genotype"]
    print(g.extend(data.columns.tolist()))
    attached_data = attached_data[["genotype"].extend(data.columns.tolist())]

    attached_genotypes = attached_data["genotype"].unique().to_list()
    if all(map(lambda s: s in DEFAULT_GENOTYPES, attached_genotypes)):
        # make genotypes be HOM, HET, WT if that seems right
        # this sucks so bad in python
        df = (
            attached_data.set_index("genotype")
            .reindex(pd.Series(DEFAULT_GENOTYPES))
            .reset_index()
        )
    else:
        # otherwise leave them like they are
        df = attached_data

    return df.sort_values("genotype")
