#!/usr/bin/env uv run
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "matplotlib",
#     "numpy",
#     "opencv-python-headless",
#     "polars",
#     "scipy",
# ]
# ///
import math
import os

import cv2
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy import ndimage

CARDINAL_DIRECTIONS = ("N", "NE", "E", "SE", "S", "SW", "W", "NW")


def load_arena_map(arena_map_file: str) -> NDArray:
    arena_map = cv2.imread(arena_map_file, cv2.IMREAD_COLOR_RGB)
    cleaned_arena_map = np.where(arena_map < 8, 0, arena_map)
    return cleaned_arena_map


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


def bin_data(all_arena_positions: dict[int, NDArray], scale: float = 3) -> NDArray:
    # plate is 127.76 mm x 85.4 mm
    bins = np.zeros((math.ceil(85.4 * scale), math.ceil(127.76 * scale)))
    for arena_positions in all_arena_positions.values():
        index_like_positions = np.floor(arena_positions * scale).astype(np.intp)
        x = index_like_positions[:, 0]
        y = index_like_positions[:, 1]
        bins[y, x] += 1
    bins = ndimage.gaussian_filter(bins, scale)
    bins = np.clip(bins, a_min=0.25, a_max=None)
    bins -= np.min(bins)
    return bins


def points_in_bin(points: NDArray, bins: NDArray, scale: float = 3) -> NDArray:
    in_bin = np.zeros(points.shape[0]).astype(bool)
    for i, (x, y) in enumerate(points):
        in_bin[i] = binned_data[math.floor(y * scale), math.floor(x * scale)] > 0
    return in_bin


def cardinals(data: NDArray) -> dict[str, NDArray]:
    # fmt: off
    cardinal_data = {
        "N":  data[np.argmin(data[:, 1])],
        "NE": data[np.argmin(data[:, 1] - data[:, 0])],
        "E":  data[np.argmin(data[:, 0])],
        "SE": data[np.argmax(data[:, 1] + data[:, 0])],
        "S":  data[np.argmax(data[:, 1])],
        "SW": data[np.argmax(data[:, 1] - data[:, 0])],
        "W":  data[np.argmax(data[:, 0])],
        "NW": data[np.argmin(data[:, 1] + data[:, 0])],
    }
    # fmt: on
    return cardinal_data


# fmt: off
assay_to_arena_type = {
    "distance_traveled":          "48_well",
    "light_dark_preference_3wpf": "12_well",
    "light_dark_preference_6dpf": "12_well",
    "light_dark_transition":      "48_well",
    "mirror_biting":              "mirror_biting",
    "sleep":                      "48_well",
    "social_preference":          "social_preference",
    "startle_response":           "48_well",
    "ymaze_4":                    "ymaze_4",
    "ymaze_15":                   "ymaze_15",
}
# fmt: on
arena_to_assay_type = {
    arena_type: set() for arena_type in set(assay_to_arena_type.values())
}
for assay_type, arena_type in assay_to_arena_type.items():
    arena_to_assay_type[arena_type].add(assay_type)
# fmt: off
arena_type_to_num_arenas = {
    "48_well":           48,
    "12_well":           12,
    "mirror_biting":     20,
    "social_preference": 10,
    "ymaze_4":            4,
    "ymaze_15":          15,
}
# some arena types end up with degenerate coordinates in certain
# directions. we should ignore those directions when matching points
arena_type_to_bad_directions = {
    "48_well":           set(),
    "12_well":           set(),
    "mirror_biting":     {"N", "E", "S", "W"},
    "social_preference": {"N", "E", "S", "W"},
    "ymaze_4":           {"N", "S"},
    "ymaze_15":          {"N", "S"},
}
# fmt: on

# create arena maps and camera parameters
for arena_type, assay_types in arena_to_assay_type.items():
    print("\n" + arena_type.upper())
    # xy schema for all xy data in this arena type
    column_names = [
        f"{direction}_A{i + 1}"
        for i in range(arena_type_to_num_arenas[arena_type])
        for direction in ("X", "Y")
    ]
    xy_schema = pl.Schema({col_name: pl.Float64 for col_name in column_names})

    # gather dataframes
    all_xy_data = []
    for assay_type in assay_types:
        directory = os.path.join("data", "all_" + assay_type)
        if os.listdir(directory):
            all_xy_data.append(
                pl.read_csv(
                    os.path.join(directory, "*"),
                    columns=column_names,
                    schema_overrides=xy_schema,
                )
            )
    if not all_xy_data:
        continue

    # combine into one long dataframe
    xy_df = pl.DataFrame(schema=xy_schema)
    for df in all_xy_data:
        xy_df.vstack(df, in_place=True)

    # get points for each arena separately
    arena_points = {}
    for i in range(arena_type_to_num_arenas[arena_type]):
        arena_idx = i + 1
        arena_points[arena_idx] = (
            xy_df.select([f"X_A{arena_idx}", f"Y_A{arena_idx}"])
            .drop_nulls()
            .unique()
            .to_numpy()
        )

    # load zantiks assay map asset
    asset_map = load_arena_map(
        os.path.join("data", "arenas", arena_type, arena_type + ".bmp")
    )

    # make binned data map to throw out tracking errors and unused arenas
    binned_data = bin_data(arena_points)
    used_arenas_indices = []
    used_arenas = np.zeros((asset_map.shape[0], asset_map.shape[1]))
    for i in range(arena_type_to_num_arenas[arena_type]):
        if len(arena_points[i + 1]) > 1000: # need > 1000 points to make a good map
            used_arenas_indices.append(i + 1)
            used_arenas += get_arena_coords(asset_map, i + 1)
    plt.figure(dpi=175)
    plt.imshow(binned_data, cmap="magma")
    plt.axis("off")
    plt.title(f"Binned Data for {arena_type}")
    plt.show()

    # get cardinal points for each arena
    good_directions = [
        direction
        for direction in CARDINAL_DIRECTIONS
        if direction not in arena_type_to_bad_directions[arena_type]
    ]
    measured_cardinals = np.zeros(
        (len(used_arenas_indices) * len(good_directions), 2),
        dtype=np.float32,
    )
    asset_cardinals = np.zeros(
        (len(used_arenas_indices) * len(good_directions), 2),
        dtype=np.float32,
    )
    for i, arena_idx in enumerate(used_arenas_indices):
        # first, all the points measured in the actual assays
        measured_points = arena_points[arena_idx]
        valid_points = measured_points[points_in_bin(measured_points, binned_data)]
        all_measured_cardinals = cardinals(valid_points)
        measured_cardinals[
            i * len(good_directions) : (i + 1) * len(good_directions)
        ] = [all_measured_cardinals[direction] for direction in good_directions]

        # then, the points from the zantiks arena map asset
        asset_points = np.flip(np.argwhere(get_arena_coords(asset_map, arena_idx)))
        all_asset_cardinals = cardinals(asset_points)
        asset_cardinals[i * len(good_directions) : (i + 1) * len(good_directions)] = [
            all_asset_cardinals[direction] for direction in good_directions
        ]

    # finally, FIND THE CAMERA PARAMETERS
    obj_points = [
        np.column_stack(
            (measured_cardinals, np.zeros(measured_cardinals.shape[0]))
        ).astype(np.float32)
    ]
    image_points = [asset_cardinals.astype(np.float32)]
    ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(
        obj_points, image_points, asset_map.shape[:2][::-1], None, None
    )

    # SAVE THEM!!!!!!
    np.savez_compressed(
        os.path.join("data", "arenas", arena_type, arena_type),
        allow_pickle=False,
        mtx=mtx,
        dist=dist,
        rvec=rvecs[0],
        tvec=tvecs[0],
    )

    # SEE THEM!!!!
    print(f"{ret=}")
    print(f"{mtx=}")
    print(f"{dist=}")
    print(f"{rvecs=}")
    print(f"{tvecs=}")

    out_points, _ = cv2.projectPoints(obj_points[0], rvecs[0], tvecs[0], mtx, dist)
    projected_points = np.squeeze(out_points)
    plt.figure(dpi=175)
    plt.scatter(
        asset_cardinals[:, 0],
        asset_cardinals[:, 1],
        alpha=0.75,
        label="Points from Asset Maps",
    )
    plt.scatter(
        projected_points[:, 0],
        projected_points[:, 1],
        alpha=0.75,
        label="Measured Points Reprojection",
    )
    plt.axis("off")
    plt.imshow(used_arenas, cmap="binary")
    plt.title(f"Reprojection of {arena_type}")
    plt.legend()
    plt.show()
