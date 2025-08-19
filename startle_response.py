import marimo

__generated_with = "0.11.22"
app = marimo.App(width="medium", css_file="data/marimo.css")


@app.cell
def _():
    import marimo as mo
    import polars as pl

    import data_utils
    return data_utils, mo, pl


@app.cell
def _(data_utils):
    xy_startle_config = data_utils.ConfigParser().parse("configs/jip3_test/6dpf/startle_response.toml")
    xy_startle_files = [
        data_utils.ZantiksFile(xy_startle_config.data.files_prefix + filepath)
        for filepath in xy_startle_config.data.files
    ]
    xy_startle_data = [data_utils.ZantiksData(startle_file) for startle_file in xy_startle_files]
    return xy_startle_config, xy_startle_data, xy_startle_files


@app.cell
def _(xy_startle_data):
    xy_startle_data[0].data
    return


@app.cell
def _(data_utils):
    startle_config = data_utils.ConfigParser().parse(
        "configs/jip3_test/6dpf/startle_response_dual.toml"
    )
    startle_files = [
        data_utils.ZantiksFile(startle_config.data.files_prefix + filepath)
        for filepath in startle_config.data.files
    ]
    startle_data = [data_utils.ZantiksData(startle_file) for startle_file in startle_files]
    return startle_config, startle_data, startle_files


@app.cell
def _(startle_data):
    startle_data[0].data
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Matching Hindges & Read
        - prepulse looks "backwards"
        - remove wells with 0 mm locomotion +/- 500 ms around the startle vibration
        - $\%PPI = \frac{\text{Startle alone} −  \text{Prepulse \& Startle}}{\text{Startle alone}} * 100$
        -  “Startle alone” was the mean locomotion at +100 ms in startle alone trials and “Prepulse/Startle” is the mean locomotion at +100 ms in prepulse/startle trials
        """
    )
    return


@app.cell
def _(pl):
    def prepare_for_prism(xy_dfs: list[pl.DataFrame], dfs: list[pl.DataFrame]):
        # cut down xy data to 1 second neighborhood around startle/prepulse timepoints
        dfs_interval_only = []
        for i, df in enumerate(dfs):
            # make intervals and cut down xy dataframes
            intervaled_dfs = []
            for startle_type in ("STARTLE", "PREPULSE"):
                # get times when startle and prepulse happen
                # subtract 1 second since the data is recorded 1 second after the startle happens
                times = (
                    (df.filter(pl.col("PHASE") == startle_type).select("RUNTIME").to_series() - 1)
                    .unique()
                    .sort()
                    .to_list()
                )

                intervals_data = []
                for j, time_point in enumerate(times):
                    intervals_data.append({"time_point": time_point, "phase_id": f"{startle_type}_{j}"})
                intervals_df = pl.DataFrame(intervals_data)
                filtered_df = (
                    xy_dfs[i]
                    .join(intervals_df, how="cross")
                    .filter(
                        (pl.col("RUNTIME") >= pl.col("time_point") - 1.1)
                        & (pl.col("RUNTIME") <= pl.col("time_point") + 1.2)
                    )
                )
                # also make 100ms time bins
                df_with_distance = (
                    filtered_df.with_columns(
                        [
                            # Create relative time from start of interval
                            (pl.col("RUNTIME") - pl.col("time_point")).alias("relative_time"),
                            # Create 100ms bins
                            (
                                ((pl.col("RUNTIME") - pl.col("time_point")) * 10)
                                .round()
                                .cast(pl.Int32)
                                * 100
                            ).alias("bin_100ms"),
                        ]
                    )
                    .sort(["ARENA", "phase_id", "RUNTIME"])
                    .with_columns(
                        [
                            # Calculate distance from previous point
                            (
                                (
                                    (pl.col("X") - pl.col("X").shift(1)) ** 2
                                    + (pl.col("Y") - pl.col("Y").shift(1)) ** 2
                                )
                                ** 0.5
                            )
                            .fill_null(0)
                            .alias("distance_step")
                        ]
                    )
                )
                df_with_distance = (
                    df_with_distance.group_by(["ARENA", "phase_id", "bin_100ms", "Cluster"])
                    .agg(
                        [
                            pl.sum("distance_step").alias("distance_100ms_total"),
                            (pl.max("RUNTIME") - pl.min("RUNTIME")).alias("actual_duration"),
                            pl.first("time_point").alias("time_point"),
                            pl.mean("X").alias("X_mean"),
                            pl.mean("Y").alias("Y_mean"),
                            pl.count().alias("n_points"),
                        ]
                    )
                    .with_columns(
                        [
                            # Normalize by actual time duration to get velocity
                            (pl.col("distance_100ms_total") / pl.col("actual_duration")).alias(
                                "velocity_100ms"
                            ),
                            # Scale to standard 100ms for comparison
                            (pl.col("distance_100ms_total") * 0.1 / pl.col("actual_duration")).alias(
                                "distance_100ms_normalized"
                            ),
                        ]
                    )
                    .sort(["ARENA", "phase_id", "bin_100ms"])
                    .drop("time_point")
                )

                intervaled_dfs.append(df_with_distance)
            dfs_interval_only.append(pl.concat(intervaled_dfs))

        # calculate how many of each genotype and add experiment ID
        wildtypes, hets, homs = 0, 0, 0
        for i in range(len(dfs_interval_only)):
            dfs_interval_only[i] = dfs_interval_only[i].with_columns(pl.lit(i).alias("Experiment"))
            genotypes = dfs[i].group_by("ARENA").agg(pl.col("Cluster").get(0))
            wildtypes += sum(genotypes["Cluster"] == "WT")
            hets += sum(genotypes["Cluster"] == "HET")
            homs += sum(genotypes["Cluster"] == "HOM")
        num_genotype_cols = max(wildtypes, hets, homs)
        print(f"You will need {num_genotype_cols} replicants in Prism!")

        # clean up to distance traveled per time point
        final_dfs = []
        for df in dfs_interval_only:
            df = df.with_columns(pl.col("phase_id").str.split("_").list.get(0).alias("startle_type"))
            df = df.group_by(["ARENA", "startle_type", "bin_100ms", "Cluster", "Experiment"]).agg(
                pl.mean("distance_100ms_total")
            )
            df = df.filter(pl.col("Cluster").is_in(("WT", "HET", "HOM")))
            final_dfs.append(df)

        combined_unfiltered = pl.concat(final_dfs).with_columns(
            (
                pl.col("Cluster")
                + "_"
                + pl.col("ARENA").cast(str)
                + "_"
                + pl.col("Experiment").cast(str)
            ).alias("Replicant ID")
        )

        results = []
        for startle_type in ("STARTLE", "PREPULSE"):
            combined = combined_unfiltered.filter(pl.col("startle_type") == startle_type)
            pivoted = combined.pivot("Replicant ID", index="bin_100ms", values="distance_100ms_total")

            def get_genotype_columns(geno: str):
                return [col for col in pivoted.columns if col.startswith(geno + "_")]

            def pad_genotype_columns(cols: list[str]):
                padded = [pivoted[col] for col in cols]
                n_pad = num_genotype_cols - len(cols)
                if n_pad > 0:
                    padded += [
                        pl.Series(
                            name=f"{cols[0].split('_')[0]}_pad{i}", values=[None] * pivoted.height
                        )
                        for i in range(n_pad)
                    ]
                return padded

            wt_cols = pad_genotype_columns(get_genotype_columns("WT"))
            het_cols = pad_genotype_columns(get_genotype_columns("HET"))
            hom_cols = pad_genotype_columns(get_genotype_columns("HOM"))

            results.append(
                pl.DataFrame([pivoted["bin_100ms"], *wt_cols, *het_cols, *hom_cols])
                .sort("bin_100ms")
                .filter((pl.col("bin_100ms") < 1000) & (pl.col("bin_100ms") > -1000))
            )

        return results
    return (prepare_for_prism,)


@app.cell
def _(prepare_for_prism, startle_data, xy_startle_data):
    startle, prepulse = prepare_for_prism(
        [data.data for data in xy_startle_data], [data.data for data in startle_data]
    )
    startle.write_csv("data/output/startle.csv")
    prepulse.write_csv("data/output/prepulse.csv")
    startle
    return prepulse, startle


if __name__ == "__main__":
    app.run()
