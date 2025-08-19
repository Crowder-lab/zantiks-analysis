#!/usr/bin/env uv run

import glob
import io
import re
import tomllib
from abc import ABC, abstractmethod
from typing import Iterable, Literal, Self, Sequence

import numpy as np
import polars as pl

type VariableDate = Sequence[int | None]


class Assay(ABC):
    @property
    @abstractmethod
    def name(self) -> str:
        pass

    @property
    @abstractmethod
    def pretty_name(self) -> str:
        pass

    @property
    @abstractmethod
    def age(self) -> str:
        pass

    @property
    def config_path(self) -> str:
        return f"configs/jip3_test/{self.age}/{self.name}.toml"

    @property
    @abstractmethod
    def total_arenas(self) -> int:
        pass

    @property
    @abstractmethod
    def arenas_per_group(self) -> int | None:
        pass

    @property
    def arena_bmp_path(self) -> str:
        return f"data/arenas/{self.name}_arenas.bmp"

    @property
    def camera_parameters_path(self) -> str:
        return f"data/camera_params/{self.name}.npz"

    @abstractmethod
    def make_tabular(self, df: pl.DataFrame) -> pl.DataFrame:
        pass


class LightDarkPreference3wpf(Assay):
    @property
    def name(self) -> str:
        return "light_dark_preference_3wpf"

    @property
    def pretty_name(self) -> str:
        return "Light/Dark Preference 3wpf"

    @property
    def age(self) -> str:
        return "3wpf"

    @property
    def total_arenas(self) -> int:
        return 12

    @property
    def arenas_per_group(self) -> int | None:
        return 4


class LightDarkPreference6dpf(Assay):
    @property
    def name(self) -> str:
        return "light_dark_preference_6dpf"

    @property
    def pretty_name(self) -> str:
        return "Light/Dark Preference 6dpf"

    @property
    def age(self) -> str:
        return "6dpf"

    @property
    def total_arenas(self) -> int:
        return 12

    @property
    def arenas_per_group(self) -> int | None:
        return 12


class LightDarkTransition(Assay):
    @property
    def name(self) -> str:
        return "light_dark_transition"

    @property
    def pretty_name(self) -> str:
        return "Light/Dark Transition"

    @property
    def age(self) -> str:
        return "6dpf"

    @property
    def total_arenas(self) -> int:
        return 48

    @property
    def arenas_per_group(self) -> int | None:
        return None


class MirrorBiting(Assay):
    @property
    def name(self) -> str:
        return "mirror_biting"

    @property
    def pretty_name(self) -> str:
        return "Mirror Biting"

    @property
    def age(self) -> str:
        return "3wpf"

    @property
    def total_arenas(self) -> int:
        return 20

    @property
    def arenas_per_group(self) -> int | None:
        return 4


class Sleep(Assay):
    @property
    def name(self) -> str:
        return "sleep"

    @property
    def pretty_name(self) -> str:
        return "Sleep"

    @property
    def age(self) -> str:
        return "6dpf"

    @property
    def total_arenas(self) -> int:
        return 48

    @property
    def arenas_per_group(self) -> int | None:
        return None

    def make_tabular(self, df: pl.DataFrame) -> pl.DataFrame:
        pattern = re.compile(r"A(\d+)_Z(\d+)")
        df = df.unpivot(
            index=("TIME", "CONDITION", "BIN_NUM"),
            on=[col for col in df.columns if re.match(pattern, col)],
            variable_name="LOCATION_CODE",
            value_name="DISTANCE",
        )
        df = df.with_columns(
            [
                pl.col("LOCATION_CODE")
                .str.extract(r"A(\d+)", 1)
                .cast(pl.Int32)
                .alias("ARENA"),
                pl.col("LOCATION_CODE")
                .str.extract(r"Z(\d+)", 1)
                .cast(pl.Int32)
                .alias("ZONE"),
            ]
        )
        df = df.drop("LOCATION_CODE")

        zero_arenas = (
            df.group_by("ARENA")
            .agg(pl.sum("DISTANCE").alias("TOTAL_DISTANCE"))
            .filter(pl.col("TOTAL_DISTANCE") == 0)
            .select("ARENA")
        )
        if not zero_arenas.is_empty():
            zero_arena_list = zero_arenas["ARENA"].to_list()
            result = df.filter(~pl.col("ARENA").is_in(zero_arena_list))
        else:
            result = df

        result = result.select(
            ("TIME", "BIN_NUM", "CONDITION", "ARENA", "ZONE", "DISTANCE")
        )
        return result


class SocialPreference(Assay):
    @property
    def name(self) -> str:
        return "social_preference"

    @property
    def pretty_name(self) -> str:
        return "Social Preference"

    @property
    def age(self) -> str:
        return "3wpf"

    @property
    def total_arenas(self) -> int:
        return 10

    @property
    def arenas_per_group(self) -> int | None:
        return 4


class StartleResponse(Assay):
    @property
    def name(self) -> str:
        return "startle_response"

    @property
    def pretty_name(self) -> str:
        return "Startle Response + Pre-pulse Inhibition"

    @property
    def age(self) -> str:
        return "6dpf"

    @property
    def total_arenas(self) -> int:
        return 48

    @property
    def arenas_per_group(self) -> int | None:
        return None

    def make_tabular(self, df: pl.DataFrame) -> pl.DataFrame:
        pattern = re.compile(r"A(\d+)")
        df = df.unpivot(
            index=["RUNTIME", "PHASE"],
            on=[col for col in df.columns if re.match(pattern, col)],
            variable_name="ARENA_ID",
            value_name="DISTANCE",
        )
        df = df.with_columns(
            [
                pl.col("ARENA_ID")
                .str.extract(r"A(\d+)", 1)
                .cast(pl.Int32)
                .alias("ARENA"),
            ]
        )
        df = df.drop("ARENA_ID")

        zero_arenas = (
            df.group_by("ARENA")
            .agg(pl.sum("DISTANCE").alias("TOTAL_DISTANCE"))
            .filter(pl.col("TOTAL_DISTANCE") == 0)
            .select("ARENA")
        )
        if not zero_arenas.is_empty():
            zero_arena_list = zero_arenas["ARENA"].to_list()
            result = df.filter(~pl.col("ARENA").is_in(zero_arena_list))
        else:
            result = df

        result = result.select(("RUNTIME", "PHASE", "ARENA", "DISTANCE"))
        return result


class StartleResponseXy(Assay):
    @property
    def name(self) -> str:
        return "startle_response"

    @property
    def pretty_name(self) -> str:
        return "Startle Response + Pre-pulse Inhibition"

    @property
    def age(self) -> str:
        return "6dpf"

    @property
    def total_arenas(self) -> int:
        return 48

    @property
    def arenas_per_group(self) -> int | None:
        return None

    def make_tabular(self, df: pl.DataFrame) -> pl.DataFrame:
        pattern = re.compile(r"([XY])_A(\d+)")
        df = df.unpivot(
            index="RUNTIME",
            on=[col for col in df.columns if re.match(pattern, col)],
            variable_name="AXIS_ARENA_ID",
            value_name="POSITION",
        )
        df = df.with_columns(
            [
                pl.col("AXIS_ARENA_ID").str.extract(r"([XY])", 1).alias("AXIS"),
                pl.col("AXIS_ARENA_ID")
                .str.extract(r"A(\d+)", 1)
                .cast(pl.Int32)
                .alias("ARENA"),
            ]
        )
        df = df.drop("AXIS_ARENA_ID")

        zero_arenas = (
            df.group_by("ARENA")
            .agg(pl.sum("POSITION").alias("TOTAL_POSITION"))
            .filter(pl.col("TOTAL_POSITION") == 0)
            .select("ARENA")
        )
        if not zero_arenas.is_empty():
            zero_arena_list = zero_arenas["ARENA"].to_list()
            result = df.filter(~pl.col("ARENA").is_in(zero_arena_list))
        else:
            result = df

        result = result.select(("RUNTIME", "ARENA", "AXIS", "POSITION"))
        return result


class Ymaze15(Assay):
    @property
    def name(self) -> str:
        return "ymaze_15"

    @property
    def pretty_name(self) -> str:
        return "Y-maze 15"

    @property
    def age(self) -> str:
        return "6dpf"

    @property
    def total_arenas(self) -> int:
        return 15

    @property
    def arenas_per_group(self) -> int | None:
        return 12


class Ymaze4(Assay):
    @property
    def name(self) -> str:
        return "ymaze_4"

    @property
    def pretty_name(self) -> str:
        return "Y-maze 4"

    @property
    def age(self) -> str:
        return "3wpf"

    @property
    def total_arenas(self) -> int:
        return 4

    @property
    def arenas_per_group(self) -> int | None:
        return 4


class ConfigParser(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, d: dict = {}):
        for key, value in d.items():
            if hasattr(value, "keys"):
                value = ConfigParser(value)
            self[key] = value

    def parse(self, path):
        with open(path, "rb") as f:
            self = ConfigParser(tomllib.load(f))

        return self


class ZantiksFile:
    assay_types = {
        "light_dark_preference_3wpf": LightDarkPreference3wpf,
        "light_dark_preference_6dpf": LightDarkPreference6dpf,
        "light_dark_transition": LightDarkTransition,
        "mirror_biting": MirrorBiting,
        "sleep": Sleep,
        "social_preference": SocialPreference,
        "startle_response": StartleResponse,
        "ymaze_15": Ymaze15,
        "ymaze_4": Ymaze4,
    }

    def _startle_response_picker(self):
        if self.is_xy:
            return StartleResponseXy
        else:
            return StartleResponse

    def __init__(self, path: str):
        self.path = path
        parts = path.split("/")
        self.filename = parts[-1]
        filename_parts = self.filename.split("-")
        self.basename = (
            parts[-2] + "/" + filename_parts[0] + "-" + filename_parts[1][:12]
        )
        self.groups = self._parse_groups(parts[-2])
        if not all(char in {"A", "B", "C", "D"} for char in self.groups):
            self.groups = None
        self.is_xy = "xy" in self.filename or "XY" in self.filename
        self.assay_type: Assay = self.assay_types[filename_parts[0]]()
        self.year = filename_parts[1][:4]
        self.month = filename_parts[1][4:6]
        self.day = filename_parts[1][6:8]

        self._parse_config()

    def __repr__(self):
        return f"ZantiksFile({self.year}, {self.month}, {self.day}, {self.assay_type}, {self.groups}, {self.filename}) - is_xy: {self.is_xy}"

    def _parse_groups(self, groups: str) -> tuple[str, ...]:
        if groups.startswith("zantiks"):
            return (groups[-1].upper(),)
        else:
            return tuple(groups.upper())

    def _parse_config(self):
        with open(self.assay_type.config_path, "rb") as f:
            config = tomllib.load(f)

        file_idx = tuple(
            map(lambda s: self.basename in s, config["data"]["files"])
        ).index(True)
        self.genotypes_path = (
            config["genotypes"]["files_prefix"] + config["genotypes"]["files"][file_idx]
        )
        self.counting_direction = config["genotypes"]["counting_directions"][file_idx]
        self.fish_used_path = (
            config["fish_used"]["files_prefix"] + config["fish_used"]["files"][file_idx]
        )


class ZantiksData:
    def __init__(self, zantiks_file: ZantiksFile, use_genotypes: bool = True):
        self.info = zantiks_file

        with open(zantiks_file.path, "r") as f:
            if not zantiks_file.is_xy:
                lines = f.readlines()
                f = io.StringIO("".join(lines[3:-1]))
            data = pl.read_csv(f)
        if zantiks_file.is_xy:
            self.data = self._expand_arenas(data)
        else:
            self.data = zantiks_file.assay_type.make_tabular(data)

        if use_genotypes:
            self._attach_genotypes()
            self.genotypes = True
        else:
            self.genotypes = None

    def __repr__(self):
        return "\n" + str(self.genotypes) + "\n" + str(self.data)

    def _expand_zones_and_arenas(self, df: pl.DataFrame) -> pl.DataFrame:
        pattern = re.compile(r"T\.A(\d+)\.Z(\d+)")

        long_df = df.unpivot(
            index=["TIME", "BIN_NUMBER"],
            on=[col for col in df.columns if re.match(pattern, col)],
            variable_name="ZONE_CODE",
            value_name="TIME_IN_ZONE",
        )
        long_df = long_df.with_columns(
            [
                pl.col("ZONE_CODE")
                .str.extract(r"T\.A(\d+)", 1)
                .cast(pl.Int32)
                .alias("ARENA"),
                pl.col("ZONE_CODE")
                .str.extract(r"Z(\d+)", 1)
                .cast(pl.Int32)
                .alias("ZONE"),
            ]
        )
        long_df = long_df.drop("ZONE_CODE")

        zero_arenas = (
            long_df.group_by("ARENA")
            .agg(pl.sum("TIME_IN_ZONE").alias("TOTAL_TIME"))
            .filter(pl.col("TOTAL_TIME") == 0)
            .select("ARENA")
        )
        if not zero_arenas.is_empty():
            zero_arena_list = zero_arenas["ARENA"].to_list()
            result = long_df.filter(~pl.col("ARENA").is_in(zero_arena_list))
        else:
            result = long_df

        result = result.select(["TIME", "BIN_NUMBER", "ARENA", "ZONE", "TIME_IN_ZONE"])
        return result

    def _expand_arenas(self, df: pl.DataFrame) -> pl.DataFrame:
        pattern = re.compile(r"[XY]_A(\d+)")

        long_df = df.unpivot(
            index="RUNTIME", on=[col for col in df.columns if re.match(pattern, col)]
        )
        long_df = long_df.with_columns(
            [
                pl.col("variable")
                .str.extract(r"(\d+)", 1)
                .cast(pl.Int32)
                .alias("ARENA"),
                pl.col("variable").str.extract(r"([XY])", 1).alias("direction"),
            ]
        )
        long_df = long_df.filter(pl.col("direction") == "X").join(
            long_df.filter(pl.col("direction") == "Y"),
            on=("ARENA", "RUNTIME"),
        )
        long_df = long_df.select(
            "RUNTIME",
            "ARENA",
            X=pl.col("value").cast(pl.Float64),
            Y=pl.col("value_right").cast(pl.Float64),
        )

        zero_arenas = (
            long_df.group_by("ARENA")
            .sum()
            .filter(pl.col("X") <= 0.0)
            .filter(pl.col("Y") <= 0.0)
            .select("ARENA")
        )
        if not zero_arenas.is_empty():
            zero_arena_list = zero_arenas["ARENA"].to_list()
            result = long_df.filter(~pl.col("ARENA").is_in(zero_arena_list))
        else:
            result = long_df

        return result

    def _attach_genotypes(self) -> None:
        # fish used data
        with open(self.info.fish_used_path, "r") as f:
            fish_used_data = f.readlines()
        fish_used_data = np.asarray([line.split() for line in fish_used_data])
        label_matrix = np.asarray(
            [
                "".join((letter, number))
                for letter in ("A", "B", "C", "D", "E", "F", "G", "H")
                for number in (
                    "01",
                    "02",
                    "03",
                    "04",
                    "05",
                    "06",
                    "07",
                    "08",
                    "09",
                    "10",
                    "11",
                    "12",
                )
            ]
        ).reshape((8, 12))
        wells_used = label_matrix[(fish_used_data == "x") | (fish_used_data == "X")]

        # genotype data
        with open(self.info.genotypes_path, "r") as f:
            genotype_data = pl.read_csv(f, columns=("Well", "Cluster"))
        genotype_data = genotype_data.with_columns(
            pl.col("Well").str.extract(r"([A-H])").alias("row"),
            pl.col("Well").str.extract(r"([0-9]+)").cast(int).alias("column"),
        )
        genotype_data = genotype_data.filter(pl.col("Well").is_in(wells_used))
        if self.info.counting_direction == "across":
            genotype_data = genotype_data.sort("row", "column")
        else:
            genotype_data = genotype_data.sort("column", "row")

        genotype_data = genotype_data.with_row_index(offset=1).select(
            "index", "Cluster"
        )
        self.data = self.data.join(
            genotype_data, left_on="ARENA", right_on="index", maintain_order="left"
        )

    def get_genotype(self, genotype: str | Iterable[str]) -> pl.DataFrame:
        if not self.genotypes:
            self._attach_genotypes()

        if isinstance(genotype, str):
            genotype = (genotype,)

        return self.data.filter(pl.col("Cluster").is_in(genotype))


class DataLoader:
    def __init__(self):
        self.pathnames: list[str] = []
        self.zantiks_files = []

    def add_by_filter(
        self,
        *,
        dates: Iterable[VariableDate] | None = None,
        assay_types: Iterable[str] | None = None,
        groups: Iterable[str] | None = None,
        data_type: Literal["zantiks"] | Literal["position"] | None = None,
    ) -> Self:
        pathglob = "**/"

        if groups:
            pathglob += f"{self._iterable_glob(groups)}/"

        if assay_types:
            pathglob += self._iterable_glob(assay_types)
        elif pathglob[-1] != "*":
            pathglob += "*"

        if dates:
            date_patterns = []
            for date in dates:
                year = str(date[0])
                month = f"{date[1]:02d}" if date[1] is not None else ""
                day = f"{date[2]:02d}" if date[2] is not None else ""
                date_patterns.append(f"{year}{month}{day}")

            pathglob += f"-{self._iterable_glob(date_patterns)}*"
        elif pathglob[-1] != "*":
            pathglob += "*"

        if pathglob[-1] != "*":
            pathglob += "*"
        pathglob += ".csv"

        matching_paths = glob.glob(pathglob, recursive=True)
        if data_type == "zantiks":
            matching_paths = [path for path in matching_paths if "xy" not in path]
        elif data_type == "position":
            matching_paths = [path for path in matching_paths if "xy" in path]

        self.add_by_name(matching_paths)

        return self

    def _iterable_glob(self, strings: Iterable[str]) -> str:
        if len(tuple(strings)) == 1:
            return strings[0]
        else:
            return f"{{{','.join(strings)}}}"

    def add_by_glob(self, pathglob: str) -> Self:
        self.add_by_name(glob.glob(pathglob))
        return self

    def add_by_name(self, paths: str | Iterable[str]) -> Self:
        if isinstance(paths, str):
            paths = (paths,)

        self.pathnames += paths
        self.zantiks_files += [ZantiksFile(path) for path in paths]

        return self

    def load_all(self, *args, **kwargs) -> list[ZantiksData]:
        data = []
        for zantiks_file in self.zantiks_files:
            data.append(ZantiksData(zantiks_file, *args, **kwargs))

        return data


if __name__ == "__main__":
    # dl = DataLoader().add_by_filter(
    #     assay_types=("startle_response",), data_type="position"
    # )
    # dfs = dl.load_all()
    # for df in dfs:
    #     print(df.get_genotype(("WT", "HOM")))

    c = ConfigParser().parse("./configs/jip3_test/6dpf/light_dark_transition.toml")
    print(c.data.files)
