import marimo

__generated_with = "0.11.22"
app = marimo.App(width="medium", css_file="data/marimo.css")


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl
    import seaborn as sns

    import data_utils
    return data_utils, mo, np, pl, plt, sns


@app.cell
def _(data_utils):
    sleep_config = data_utils.ConfigParser().parse("configs/jip3_test/6dpf/sleep.toml")
    sleep_files = [
        data_utils.ZantiksFile(sleep_config.data.files_prefix + path)
        for path in sleep_config.data.files
    ]
    sleep_data = [data_utils.ZantiksData(sleep_file) for sleep_file in sleep_files]
    return sleep_config, sleep_data, sleep_files


@app.cell
def _(pl, sleep_data):
    for sd in sleep_data:
        print(sd.data.filter(pl.col("CONDITION") == "DARK"))
    return (sd,)


@app.cell
def _(np, pl, plt, sns):
    def make_plot(df: pl.DataFrame, title: str) -> None:
        df = df.with_columns(((pl.col("TIME") - df[0, "TIME"]).cast(pl.Int32) // 3600).alias("Hour"))
        df = df.group_by("ARENA", "Hour").agg(
            pl.col("CONDITION").get(0).eq("BRIGHT").alias("Light"),
            pl.sum("DISTANCE").alias("Distance (mm)"),
            pl.col("Cluster").get(0).alias("Genotype"),
        )
        df = df.filter(pl.col("Genotype").is_in(("WT", "HET", "HOM")))
        df = df.sort("Genotype", "Hour")
        plt.figure(dpi=200)
        sns.lineplot(
            df,
            x="Hour",
            y="Distance (mm)",
            hue="Genotype",
            hue_order=("WT", "HET", "HOM"),
        )
        plt.title(title)
        plt.xticks(np.arange(0, 24, 2))
        plt.show()
    return (make_plot,)


@app.cell
def _(make_plot, sleep_data):
    for sleep_datum in sleep_data:
        make_plot(sleep_datum.data, sleep_datum.info.filename)
    return (sleep_datum,)


@app.cell
def _(data_utils, np, pl, plt, sns):
    def make_combined_plot(data: list[data_utils.ZantiksData]):
        dfs = []
        for df in data:
            df = df.with_columns(
                ((pl.col("TIME") - df[0, "TIME"]).cast(pl.Int32) // 3600).alias("Hour")
            )
            df = df.group_by("ARENA", "Hour").agg(
                pl.col("CONDITION").get(0).eq("BRIGHT").alias("Light"),
                pl.sum("DISTANCE").alias("Distance (mm)"),
                pl.col("Cluster").get(0).alias("Genotype"),
            )
            df = df.filter(pl.col("Genotype").is_in(("WT", "HET", "HOM")))
            dfs.append(df)
        combined = pl.concat(dfs)
        combined = combined.sort("Genotype", "Hour")
        plt.figure(dpi=200)
        sns.lineplot(
            combined,
            x="Hour",
            y="Distance (mm)",
            hue="Genotype",
            hue_order=("WT", "HET", "HOM"),
        )
        plt.title("All Days")
        plt.xticks(np.arange(0, 24, 2))
        return plt.gca()
    return (make_combined_plot,)


@app.cell
def _(make_combined_plot, sleep_data):
    make_combined_plot([data.data for data in sleep_data])
    return


@app.cell
def _(pl):
    def prepare_for_prism(dfs: list[pl.DataFrame]):
        # calculate how many of each genotype and add experiment ID
        wildtypes, hets, homs = 0, 0, 0
        for i in range(len(dfs)):
            dfs[i] = dfs[i].with_columns(pl.lit(i).alias("Experiment"))
            genotypes = (
                dfs[i]
                .group_by("ARENA")
                .agg(pl.col("Cluster").get(0))
                .filter(pl.col("Cluster").is_in(("WT", "HET", "HOM")))
            )
            wildtypes += sum(genotypes["Cluster"] == "WT")
            hets += sum(genotypes["Cluster"] == "HET")
            homs += sum(genotypes["Cluster"] == "HOM")
        num_genotype_cols = max(wildtypes, hets, homs)
        print(f"You will need {num_genotype_cols} replicants in Prism!")

        # normal combining to get hourly distance
        final_dfs = []
        for df in dfs:
            df = df.with_columns(
                ((pl.col("TIME") - df[0, "TIME"]).cast(pl.Int32) // 3600).alias("Hour")
            )
            df = df.group_by("ARENA", "Hour").agg(
                pl.col("CONDITION").get(0).eq("BRIGHT").alias("Light"),
                pl.sum("DISTANCE").alias("Distance (mm)"),
                pl.col("Cluster").get(0).alias("Genotype"),
                pl.col("Experiment").get(0),
            )
            df = df.filter(pl.col("Genotype").is_in(("WT", "HET", "HOM")))
            final_dfs.append(df)

        combined = pl.concat(final_dfs).with_columns(
            (
                pl.col("Genotype")
                + "_"
                + pl.col("ARENA").cast(str)
                + "_"
                + pl.col("Experiment").cast(str)
            ).alias("Replicant ID")
        )
        pivoted = combined.pivot("Replicant ID", index="Hour", values="Distance (mm)")

        def get_genotype_columns(geno: str):
            return [col for col in pivoted.columns if col.startswith(geno + "_")]

        def pad_genotype_columns(cols: list[str]):
            padded = [pivoted[col] for col in cols]
            n_pad = num_genotype_cols - len(cols)
            if n_pad > 0:
                padded += [
                    pl.Series(name=f"{cols[0].split('_')[0]}_pad{i}", values=[None] * pivoted.height)
                    for i in range(n_pad)
                ]
            return padded

        wt_cols = pad_genotype_columns(get_genotype_columns("WT"))
        het_cols = pad_genotype_columns(get_genotype_columns("HET"))
        hom_cols = pad_genotype_columns(get_genotype_columns("HOM"))

        result = pl.DataFrame([pivoted["Hour"], *wt_cols, *het_cols, *hom_cols]).sort("Hour")

        return result
    return (prepare_for_prism,)


@app.cell
def _(prepare_for_prism, sleep_data):
    for_prism = prepare_for_prism([data.data for data in sleep_data])
    for_prism.write_csv("data/output/all_sleep.csv")
    for_prism
    return (for_prism,)


if __name__ == "__main__":
    app.run()
