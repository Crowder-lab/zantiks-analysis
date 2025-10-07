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
        for csv_suffix in csv_suffixes:
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


def get_arena_coords(arena_map: NDArray, arena: int) -> NDArray[bool]:
    color_cycle = np.asarray(((120, 120, 248),  # blue
                              (120, 248, 120),  # green
                              (120, 248, 248),  # cyan
                              (248, 120, 120),  # red
                              (248, 120, 248),  # pink
                              (248, 248, 120))) # yellow
    num_colors = 6
    darken_by = 8
    times_to_darken, color = divmod(arena - 1, num_colors)
    arena_color = color_cycle[color] - times_to_darken * darken_by
    return np.all(arena_map == arena_color, axis=-1)


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
