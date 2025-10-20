#!/usr/bin/env uv run
import math
import os

import cv2
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import tqdm
from scipy import ndimage

import utils

# assay variable setup
suffixes = {"genotypes": "_genotypes.csv", "fish_used": "_fish.txt"}
all_files = utils.find_heatmap_data("social_preference", suffixes)

# get files
main_files = all_files["main_files"]
wildtype_files = all_files["wildtype_files"]
arena_map = utils.load_arena_map(
    os.path.join("data", "arenas", "social_preference", "social_preference.bmp")
)
camera_parameters = utils.load_camera_params(
    os.path.join("data", "arenas", "social_preference", "social_preference.npz")
)

# see if there's usable wildtype data
wildtype_exists = len(wildtype_files) != 0
wildtype_only_names = set(wildtype_files.keys()).difference(set(main_files.keys()))
wildtype_should_be_analyzed = wildtype_exists and len(wildtype_only_names) > 0


def save_arena_plot(data, arena: int, title: str, save_to: str) -> None:
    # make mask
    arena_coords = utils.get_arena_coords(arena_map, arena)
    mask = np.ma.masked_array(data, mask=arena_coords)
    mask[np.logical_not(arena_coords)] = 0.0

    # crop to size
    row_min, row_max, col_min, col_max = utils.find_crop_coordinates(arena_coords)
    cropped_data = data[row_min:row_max, col_min:col_max]
    cropped_mask = mask[row_min:row_max, col_min:col_max]

    # plot and save
    plt.figure(dpi=100)
    plt.title(title)
    plt.imshow(cropped_data, cmap="viridis")
    plt.imshow(cropped_mask, cmap="binary")
    plt.axis("off")
    plt.savefig(save_to, dpi=500, bbox_inches="tight")
    plt.show()


def save_whole_plot(data, mask, title: str, save_to: str) -> None:
    plt.figure(dpi=100)
    plt.title(title)
    im = plt.imshow(data, cmap="viridis")
    plt.imshow(mask, cmap="binary")
    plt.colorbar(im, label="Frequency", orientation="vertical", shrink=0.805, aspect=13)
    plt.xticks([])
    plt.yticks([])
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.box(on=True)
    plt.savefig(save_to, dpi=500, bbox_inches="tight")
    plt.show()


def analyze(files):
    analyzed_data = {"run": {}, "by_genotype": {}}
    for prefix_name in files.keys():
        group = files[prefix_name]
        main_data = utils.load_xy(group["main"]).drop("RUNTIME").with_row_index()
        genotypes = utils.load_genotypes(
            group["genotypes"], group["fish_used"], "across"
        )

        long_x_df = (
            main_data.lazy()
            .rename({"index": "index_x"})
            .unpivot(
                on=pl.selectors.starts_with("X"),
                index="index_x",
                variable_name="arena_x",
                value_name="X",
            )
            .with_columns(
                arena_x=pl.col("arena_x").str.extract(r"A(\d+)$").cast(pl.UInt8)
            )
            .sort(["index_x", "arena_x"])
        )
        long_y_df = (
            main_data.lazy()
            .rename({"index": "index_y"})
            .unpivot(
                on=pl.selectors.starts_with("Y"),
                index="index_y",
                variable_name="arena_y",
                value_name="Y",
            )
            .with_columns(
                arena_y=pl.col("arena_y").str.extract(r"A(\d+)$").cast(pl.UInt8)
            )
            .sort(["index_y", "arena_y"])
        )
        long_df = (
            pl.concat([long_x_df, long_y_df], how="horizontal", parallel=True)
            .drop("index_x", "arena_x", "index_y", "arena_y")
            .drop_nulls()
            .with_columns(pl.lit(0.0).alias("Z"))
            .cast(pl.Float64)
            .collect()
        )

        # use camera parameters to unwarp points
        world_points = long_df.to_numpy()
        rotation_matrix, _ = cv2.Rodrigues(camera_parameters["rvec"])
        camera_points = (
            np.dot(rotation_matrix, world_points.T) + camera_parameters["tvec"]
        ).T
        points_2d = camera_points[:, :2] / camera_points[:, 2:]
        r2 = np.pow(points_2d[:, 0], 2) + np.pow(points_2d[:, 1], 2)
        r4 = np.pow(r2, 2)
        r6 = np.pow(r2, 3)
        k1, k2, p1, p2, k3 = camera_parameters["dist"][0]
        radial_distortion = 1 + k1 * r2 + k2 * r4 + k3 * r6
        x, y = points_2d.T
        tangential_distortion = np.asarray(
            [
                2 * p1 * x * y + p2 * (r2 + 2 * np.pow(x, 2)),
                2 * p2 * x * y + p1 * (r2 + 2 * np.pow(y, 2)),
            ]
        ).T
        points_2d = points_2d * radial_distortion.reshape(-1, 1) + tangential_distortion
        points = (
            np.dot(camera_parameters["mtx"][:2, :2], points_2d.T).T
            + camera_parameters["mtx"][:2, 2]
        )

        # create map
        valid_arenas = np.zeros((arena_map.shape[0], arena_map.shape[1])).astype(bool)
        heatmap = np.zeros_like(valid_arenas).astype(np.float32)
        valid_genotypes = genotypes.filter(
            (pl.col("genotype") != "<Excluded>") & pl.col("genotype").is_not_null()
        )
        for arena in valid_genotypes["arena"]:
            valid_arenas |= utils.get_arena_coords(arena_map, arena)
        for x, y in tqdm.tqdm(points):
            if x <= 0.0 or y <= 0.0 or x > heatmap.shape[1] or y > heatmap.shape[0]:
                continue
            heatmap[math.floor(y), math.floor(x)] += 1

        # create appropriate masks for data
        arena_mask = np.ma.masked_array(heatmap, mask=valid_arenas)
        arena_mask[np.logical_not(valid_arenas)] = 0.0
        plottable_heatmap = ndimage.gaussian_filter(heatmap.copy(), 10)

        analyzed_data["run"][prefix_name] = {
            "mask": arena_mask,
            "data": plottable_heatmap,
            "arenas": valid_genotypes["arena"].to_list(),
        }

        # now by genotype
        for genotype in valid_genotypes["genotype"].unique():
            # initialize dictionary if necessary
            if genotype not in analyzed_data["by_genotype"]:
                analyzed_data["by_genotype"][genotype] = {}

            # recalculate using only arenas in the genotype
            valid_arenas = np.zeros((arena_map.shape[0], arena_map.shape[1])).astype(
                bool
            )
            for arena in valid_genotypes.filter(pl.col("genotype") == genotype)[
                "arena"
            ]:
                valid_arenas |= utils.get_arena_coords(arena_map, arena)
            arena_mask = np.ma.masked_array(heatmap, mask=valid_arenas)

            analyzed_data["by_genotype"][genotype][prefix_name] = {
                "mask": arena_mask,
                "data": plottable_heatmap,
            }

    return analyzed_data


main_data = analyze(main_files)
for prefix, data in main_data["run"].items():
    save_whole_plot(
        data["data"].copy(),
        data["mask"],
        prefix,
        os.path.join("data", "social_preference", "output", prefix),
    )
    for arena in data["arenas"]:
        save_arena_plot(
            data["data"].copy(),
            arena,
            f"{prefix}: Arena {arena}",
            os.path.join("data", "social_preference", "output", f"{prefix}_{arena}"),
        )
for genotype, prefixes in main_data["by_genotype"].items():
    for prefix, data in prefixes.items():
        save_whole_plot(
            data["data"].copy(),
            data["mask"],
            f"{prefix}: {genotype}",
            os.path.join("data", "social_preference", "output", f"{prefix}_{genotype}"),
        )
