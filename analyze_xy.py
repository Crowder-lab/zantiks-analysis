import marimo

__generated_with = "0.18.0"
app = marimo.App(
    width="medium",
    css_file="/Users/dunb/.config/marimo/default.css",
)


@app.cell
def _():
    import itertools

    import marimo as mo
    import matplotlib.pyplot as plt
    import polars as pl
    import seaborn as sns
    return itertools, pl, plt, sns


@app.cell
def _():
    xy_file_path = "data/distance_traveled/ZC4H2_2025-11-04_A_xy.csv"
    return (xy_file_path,)


@app.cell
def _(itertools):
    zantiks_indices = tuple(
        map(
            lambda p: p[0] + format(p[1], "02d"),
            itertools.product(("A", "B", "C", "D", "E", "F"), range(1, 9)),
        )
    )
    return


@app.cell
def _(pl, xy_file_path):
    wide_df = pl.read_csv(xy_file_path).drop("RUNTIME").with_row_index()
    long_x_df = (
        wide_df.lazy()
        .rename({"index": "index_x"})
        .unpivot(
            on=pl.selectors.starts_with("X"),
            index="index_x",
            variable_name="arena_x",
            value_name="X",
        )
        .with_columns(arena_x=pl.col("arena_x").str.extract(r"A(\d+)$").cast(pl.UInt8))
        .sort(["index_x", "arena_x"])
    )
    long_y_df = (
        wide_df.lazy()
        .rename({"index": "index_y"})
        .unpivot(
            on=pl.selectors.starts_with("Y"),
            index="index_y",
            variable_name="arena_y",
            value_name="Y",
        )
        .with_columns(arena_y=pl.col("arena_y").str.extract(r"A(\d+)$").cast(pl.UInt8))
        .sort(["index_y", "arena_y"])
    )
    df = (
        pl.concat([long_x_df, long_y_df], how="horizontal", parallel=True)
        .drop("index_x", "arena_y", "index_y")
        .cast({"X": pl.Float32, "Y": pl.Float32})
        .collect()
        .rename({"arena_x": "ARENA"})
        .sort(by="ARENA")
    )
    df
    return (df,)


@app.cell
def _(df, pl):
    missing = (
        df.group_by("ARENA")
        .agg(
            ((pl.col("X").is_null() & pl.col("Y").is_null()).sum() / pl.len() * 100).alias(
                "% missing values"
            )
        )
        .sort(by="% missing values", descending=True)
    )
    missing
    return


@app.cell
def _(df, pl):
    distance = df.with_columns(
        (
            ((pl.col("X") - pl.col("X").shift()) ** 2 + (pl.col("Y") - pl.col("Y").shift()) ** 2) ** 0.5
        ).alias("distance")
    )
    distance
    return (distance,)


@app.cell
def _(distance, pl, plt, sns):
    for arena in distance["ARENA"].unique().to_list():
        arena_distance = distance.filter(pl.col("ARENA") == arena)
        print(arena)
        print("\tsum\t\t", arena_distance["distance"].sum())
        print("\tmean\t", arena_distance["distance"].mean())
        print("\tstd\t\t", arena_distance["distance"].std())
        print("\tmin\t\t", arena_distance["distance"].min())
        print("\tmax\t\t", arena_distance["distance"].max())
        plt.figure(figsize=(8, 6))
        sns.histplot(
            arena_distance,
            x="distance",
            stat="proportion",
            binwidth=0.1,
            binrange=(0, 1.5),
            log_scale=True,
        )
        plt.title(arena)
        plt.show()
    return


if __name__ == "__main__":
    app.run()
