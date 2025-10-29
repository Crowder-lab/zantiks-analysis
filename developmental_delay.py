#!/usr/bin/env uv run
import math
import os

import cv2
import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy.ndimage import gaussian_filter1d

import utils

# assay variable setup
suffixes = {"picture": ".png", "genotypes": "_genotypes.csv", "fish_used": "_fish.txt"}
all_files = utils.find_data("developmental_delay", suffixes)

# get files
main_files = all_files["main_files"]
wildtype_files = all_files["wildtype_files"]
arena_map = utils.load_arena_map(
    os.path.join("data", "arenas", "96_well_embryo", "96_well_embryo.bmp")
)

# see if there's usable wildtype data
wildtype_exists = len(wildtype_files) != 0
wildtype_only_names = set(wildtype_files.keys()).difference(set(main_files.keys()))
wildtype_should_be_analyzed = wildtype_exists and len(wildtype_only_names) > 0


# functions for image processing
def normalize(stats_from: NDArray, to_normalize: NDArray | None = None) -> NDArray:
    """
    Normalize to [0, 1] based on minimum and maximum values.
    """
    nmin = np.nanmin(stats_from)
    nmax = np.nanmax(stats_from)
    if to_normalize is None:
        to_normalize = stats_from
    normalized = (to_normalize - nmin) / (nmax - nmin)
    return normalized


def monotonically_increasing(a: float, b: float, c: float) -> bool:
    # no calc early return
    if c == 0:
        return a >= 0 and a >= -2 * b and (a != 0 or b != 0)

    radicand = 4 * b * b - 12 * a * c
    # some calc early return
    if c > 0:
        if radicand < 0:
            return True
        elif b >= 0 and b + 3 * c <= 0:
            return True

    radical = math.sqrt(radicand)
    plus_root = (-2 * b + radical) / (6 * c)
    minus_root = (-2 * b - radical) / (6 * c)
    return (c < 0 and plus_root <= 0 and minus_root >= 1) or (
        c > 0 and (plus_root <= 0 or minus_root >= 1)
    )


def discrete_entropy_vectorized(
    cm_x: float,
    cm_y: float,
    max_distance: float,
    parameters: tuple[float, float, float],
    image: NDArray,
):
    a, b, c = parameters
    num_rows, num_cols = image.shape

    # make bigger than needed
    histogram = np.zeros((512,))

    # normalized distance to CoM
    r = (
        np.sqrt(
            np.pow(np.broadcast_to(np.arange(num_cols), image.shape) - cm_x, 2)
            + np.pow(
                np.broadcast_to(np.arange(num_rows), (num_cols, num_rows)).T - cm_y, 2
            )
        )
        / max_distance
    )

    # gain function and pixel luminance value
    g = 1 + a * np.pow(r, 2) + b * np.pow(r, 4) + c * np.pow(r, 6)
    intensity = image * g

    # map luminance to bins
    bin = 255 * np.log1p(intensity) / math.log(256)
    floor_bin = np.floor(bin)
    floor_bin_vals = 1 + floor_bin - bin
    ceil_bin = np.ceil(bin)
    ceil_bin_vals = ceil_bin - bin

    # sometimes luminance exceeds 255 after being adjusted
    max_histogram_bin = math.ceil(np.max(bin)) if np.max(bin) > 255 else 255

    for i in range(math.ceil(max_histogram_bin)):
        histogram[i] += np.sum(floor_bin_vals[floor_bin == i]) + np.sum(
            ceil_bin_vals[ceil_bin == i]
        )

    # smooth out histogram
    histogram = gaussian_filter1d(histogram[:max_histogram_bin], sigma=4)
    scaled_histogram = histogram / np.sum(histogram)

    # discrete entropy
    return -np.sum(
        np.where(scaled_histogram > 0, scaled_histogram * np.log(scaled_histogram), 0)
    )


def find_parameters_vectorized(
    cm_x: float, cm_y: float, max_distance: float, image: NDArray
):
    a = b = c = 0
    delta = 2
    min_h = None

    explored = set()
    min_h = math.inf
    while delta > 1 / 256:
        initial_guess = (a, b, c)
        guess_matrix = np.broadcast_to(np.asarray((a, b, c)), (3, 3))
        delta_matrix = delta * np.identity(3)
        guesses = np.hstack(
            (guess_matrix + delta_matrix, guess_matrix - delta_matrix)
        ).reshape((-1, 3))
        for guess_arr in guesses:
            guess = tuple(guess_arr.tolist())
            if guess not in explored:
                explored.add(guess)

                if monotonically_increasing(*guess):
                    current_h = discrete_entropy_vectorized(
                        cm_x, cm_y, max_distance, guess, image
                    )

                    if current_h < min_h:
                        min_h = current_h
                        a, b, c = guess

        if initial_guess == (a, b, c):
            delta /= 2

    return a, b, c


def vignetting_correction_vectorized(image: NDArray) -> NDArray:
    """
    Adapted from https://github.com/Hiroki39/Vignetting-Correction
    """
    original_image = image.copy()
    image = cv2.transform(original_image, np.array([[0.2126, 0.7152, 0.0722]]))
    num_rows, num_cols = image.shape

    # calculate centers of mass
    cm_x = np.sum(image * np.arange(num_cols)) / np.sum(image)
    cm_y = np.sum(image.T * np.arange(num_rows)) / np.sum(image)
    distance_array = np.asarray(
        ((0, 0), (0, num_rows), (num_cols, 0), (num_cols, num_rows))
    )
    max_distance = math.sqrt(
        np.max(np.sum(np.pow(distance_array - np.asarray((cm_x, cm_y)), 2), axis=1))
    )

    a, b, c = find_parameters_vectorized(cm_x, cm_y, max_distance, image)

    r = (
        np.sqrt(
            np.pow(np.broadcast_to(np.arange(num_cols), image.shape) - cm_x, 2)
            + np.pow(
                np.broadcast_to(np.arange(num_rows), (num_cols, num_rows)).T - cm_y, 2
            )
        )
        / max_distance
    )
    g = 1 + a * np.pow(r, 2) + b * np.pow(r, 4) + c * np.pow(r, 6)
    modified = original_image * np.dstack((g, g, g))
    return modified


# analysis function
def analyze(files):
    analyzed_data = {"percent_difference": {}}
    for prefix_name in files.keys():
        group = files[prefix_name]
        main_data = cv2.imread(group["picture"], cv2.IMREAD_COLOR_BGR)
        genotypes = utils.load_genotypes(
            group["genotypes"], group["fish_used"], "across"
        )

        # arenas are in a really ridiculous order--rearrange
        genotypes = genotypes.with_columns(
            (
                pl.col("row").str.encode("hex").str.to_integer(base=16, dtype=pl.UInt8)
                - 65
            ).alias("row_as_num")
        )
        # fmt: off
        genotypes = genotypes.with_columns(
            pl.when(pl.col("row_as_num").is_in({0, 1, 2, 3}) & (pl.col("column") <= 6))
              .then(pl.col("row_as_num") * 6 + pl.col("column"))

              .when(pl.col("row_as_num").is_in({0, 1, 2, 3}) & (pl.col("column") >  6))
              .then(pl.col("row_as_num") * 6 + pl.col("column") - 6 + 24)

              .when(pl.col("row_as_num").is_in({4, 5, 6, 7}) & (pl.col("column") <= 6))
              .then((pl.col("row_as_num") - 4) * 6 + pl.col("column") + 48)

              .when(pl.col("row_as_num").is_in({4, 5, 6, 7}) & (pl.col("column") >  6))
              .then((pl.col("row_as_num") - 4) * 6 + pl.col("column") - 6 + 72)

              .alias("row_id")
        )
        # fmt: on
        genotypes = genotypes.filter(pl.col("genotype").is_in(utils.DEFAULT_GENOTYPES))

        # make half size if image is old full resolution
        if main_data.shape[0] == 1080:
            main_data = cv2.resize(main_data, None, fx=0.5, fy=0.5)

        # remove vignette
        corrected = vignetting_correction_vectorized(main_data)
        grey_float = cv2.cvtColor(corrected.astype(np.float32), cv2.COLOR_BGR2GRAY)

        # clear out timestamp
        grey_float[520:, :320] = (
            np.mean(grey_float[:520, :]) * 720 * 520
            + np.mean(grey_float[520:, 320:]) * 400 * 20
        ) / 382400
        grey_float = 1.0 - normalize(grey_float)
        # print(grey_float.min(), grey_float.max(), grey_float.mean())
        # plt.figure(dpi=175)
        # plt.title(prefix_name)
        # plt.imshow(grey_float, cmap="magma")
        # plt.show()

        # filter and transform
        # grey_float = 1 / (1 + np.exp(-(2**3) * (grey_float - 0.9)))
        # grey_float = cv2.bilateralFilter(grey_float, 3, 2, 2)
        # grey_float = normalize(grey_float)
        # plt.figure(dpi=175)
        # plt.title(prefix_name)
        # plt.imshow(grey_float, cmap="magma")
        # plt.show()

        # find average value for all arenas
        filled_arenas = genotypes["row_id"].to_list()
        filled_arena_coords = np.zeros_like(grey_float, dtype=bool)
        for filled_arena in filled_arenas:
            filled_arena_coords |= utils.get_arena_coords(arena_map, filled_arena)
        masked_filled_arenas = np.ma.masked_array(
            grey_float, mask=np.logical_not(filled_arena_coords)
        )
        filled_arena_mean = masked_filled_arenas.sum() / len(filled_arenas)

        # for each arena
        arena_mean_comparisons = {"arena": [], "percent_difference": []}
        for i in filled_arenas:
            arena_coords = utils.get_arena_coords(arena_map, i)
            masked_arena = np.ma.masked_array(
                grey_float, mask=np.logical_not(arena_coords)
            )
            arena_mean_comparisons["arena"].append(i)
            arena_mean_comparisons["percent_difference"].append(
                ((masked_arena.sum() / filled_arena_mean) - 1.0) * 100
            )
        mean_comparison_data = pl.DataFrame(arena_mean_comparisons)
        analyzed_data["percent_difference"][prefix_name] = utils.attach_genotypes(
            mean_comparison_data, genotypes
        )

        # arena_only = grey_float.copy()
        # arena_only[np.logical_not(arena_coords)] = 0
        # grey = np.round(arena_only * 255).astype(np.uint8)
        # circles = cv2.HoughCircles(grey, cv2.HOUGH_GRADIENT, 1, 1e-3, param2=16, minRadius=8, maxRadius=14)
        # print(circles)
        # if circles is not None:
        #     circles = np.uint16(np.round(circles[(circles[:, :, 2] < 12) & (circles[:, :, 2] > 9)]))
        #     print(circles)
        #     if len(circles) > 0:
        #         x, y, r = circles[0]
        #         color = cv2.cvtColor(grey, cv2.COLOR_GRAY2RGB)
        #         cv2.circle(color, (x, y), r, (0, 255, 0), 1)
        #         cv2.circle(color, (x, y), 1, (255, 0, 0), 2)
        #         plt.figure(dpi=175)
        #         plt.title(prefix_name)
        #         plt.imshow(color)
        #         plt.show()

    return analyzed_data


main_data = analyze(main_files)
for prefix_name in main_files.keys():
    prism_data = utils.column_data(
        main_data["percent_difference"][prefix_name], "percent_difference"
    )
    with open(
        os.path.join(
            "data",
            "developmental_delay",
            "output",
            prefix_name + "_PERCENT-DIFFERENCE.csv",
        ),
        "w",
    ) as f:
        prism_data.write_csv(f)

# combine data (including wildtypes if possible)
if wildtype_should_be_analyzed:
    wildtype_data = analyze(wildtype_files)
    all_names = set(wildtype_files.keys()) | set(main_files.keys())

    combined_wildtype = (
        pl.concat(
            [
                df.with_columns(id=pl.lit(name))
                for name, df in wildtype_data["percent_difference"].items()
            ],
            how="vertical",
        )
        .filter(pl.col("genotype") == "WT")
        .filter(pl.col("id").is_in(wildtype_only_names))
    )
    combined_main = pl.concat(
        [
            df.with_columns(id=pl.lit(name))
            for name, df in main_data["percent_difference"].items()
        ],
        how="vertical",
    )
    all_percent_difference_data = pl.concat(
        (combined_wildtype, combined_main), how="vertical"
    )
else:
    all_names = set(main_files.keys())

    all_percent_difference_data = pl.concat(
        [
            df.with_columns(id=pl.lit(name))
            for name, df in main_data["percent_difference"].items()
        ],
        how="vertical",
    )

# analyze and save combined data
for_prism = utils.add_numbering(all_percent_difference_data, all_names, "genotype")
prism_data = utils.column_data(for_prism, "percent_difference")
with open(
    os.path.join(
        "data",
        "developmental_delay",
        "output",
        "combined_PERCENT-DIFFERENCE.csv",
    ),
    "w",
) as f:
    prism_data.write_csv(f)
