#!/usr/bin/env python
# Copyright Mat Beale and PAM Informatics, Wellcome Sanger Institute, 2024
import sys
import pandas as pd
import numpy as np
import argparse
import logging

from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parses the output of `samtools depth` to determine coverage statistics for each region in a given BED file. Note that this script currently assumes a single reference chromosome."
    )
    parser.add_argument(
        "--samtools-depth",
        "-s",
        type=validate_samtools_depth,
        required=True,
        dest="samtools_depth",
        help="path to `samtools depth` output file",
    )
    parser.add_argument(
        "--bed",
        "-b",
        type=validate_bed_path,
        required=True,
        dest="bed",
        help=r"path to BED file (TSV) defining regions. If BED3 format (<chrom>\t<start>\t<end>), arbitrary region names will be generated. If BED4 format (<chrom>\t<start>\t<end>\t<name>), the region names will be used if present and arbitrary names will be used for any missing names.",
    )
    parser.add_argument(
        "--threshold",
        "-t",
        type=validate_coverage_threshold,
        required=True,
        dest="coverage_thresholds",
        help="comma separated list of coverage thresholds",
    )
    parser.add_argument(
        "--sample-name",
        "-n",
        type=str,
        dest="sample_name",
        help="sample name to conveniently link sample name as metadata in coverage summary",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        dest="output",
        default=".",
        help="path to output directory",
    )
    parser.set_defaults(feature=False)

    return parser.parse_args()


## ARG VALIDATION
def validate_coverage_threshold(coverage_threshold_arg: str) -> list[float]:
    valid_thresholds = set()
    coverage_thresholds = coverage_threshold_arg.split(",")
    for threshold in coverage_thresholds:
        try:
            valid_threshold = float(threshold)
        except ValueError:
            logging.error(f"Given coverage threshold '{threshold}' is not a valid.")
            sys.exit(1)
        if valid_threshold < 0:
            logging.error(f"Given coverage threshold '{threshold}' cannot be negative.")
            sys.exit(1)
        valid_thresholds.add(valid_threshold)
    return list(sorted(valid_thresholds))


def validate_bed_path(path: str) -> Path:
    bed_path = Path(path)
    if not bed_path.is_file():
        logging.error(f"Given path '{path}' is not a file.")
        sys.exit(1)
    return bed_path


def validate_samtools_depth(path: str) -> Path:
    samtools_depth_path = Path(path)
    if not samtools_depth_path.is_file():
        logging.error(f"Given samtools depth output file '{path}' is not a valid file.")
        sys.exit(1)
    return samtools_depth_path


## PARSING
def parse_depth_data(depth_file: Path) -> pd.DataFrame:
    depth_df = pd.read_csv(
        depth_file, sep="\t", index_col=0, header=None, names=["Ref", "Pos", "Cov"]
    )
    # Anything that isn't covered (is na) make 0
    depth_df[pd.isnull(depth_df["Cov"])] = 0
    return depth_df


def parse_bed_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, header=None, sep="\t")
    if len(df.columns) < 3:
        logging.error(
            f"Given BED file '{path}' has fewer than 3 columns (please ensure the file has the correct format)."
        )
        sys.exit(1)
    if len(df.columns) >= 4:
        region_names = df.iloc[:, 3]
        duplicate_region_names = region_names.dropna().duplicated()
        if duplicate_region_names.any():
            logging.error(
                f"Given BED file '{path}' has duplicate region names in the 4th column (please ensure region names are unique)."
            )
            sys.exit(1)
    if df.iloc[:, 0:3].isna().values.any():
        logging.error(
            f"Given BED file '{path}' has missing values in one of the first 3 columns (please ensure these columns do not have missing values)."
        )
        sys.exit(1)
    integer_columns = (1, 2)
    col_type_error = False
    for column in integer_columns:
        if not df.iloc[:, column].dtypes == "int64":
            col_type_error = True
            logging.error(
                f"Given BED file '{path}' has non-integer values in column {column + 1}."
            )
        if (df.iloc[:, column] < 0).any():
            logging.error(
                f"Given BED file '{path}' has negative values in column {column + 1}. Please ensure this column contains only non-negative integers."
            )
    if col_type_error:
        sys.exit(1)
    df = fill_in_missing_region_names(df)
    return df


def fill_in_missing_region_names(bed_df: pd.DataFrame) -> pd.DataFrame:
    bed_df_copy = bed_df.copy()
    start = bed_df_copy.iloc[:, 1]
    end = bed_df_copy.iloc[:, 2]
    generated_region_names = "locus_" + start.astype(str) + "_" + end.astype(str)
    if len(bed_df_copy.columns) > 3:
        region_names = bed_df_copy.iloc[:, 3]
        missing_region_name_index = bed_df_copy.iloc[:, 3].isna()
        if missing_region_name_index.any():
            replacement_names = generated_region_names[missing_region_name_index]
            bed_df_copy.iloc[missing_region_name_index, 3] = replacement_names
    elif len(bed_df_copy.columns) == 3:
        bed_df_copy[3] = generated_region_names
    else:
        raise ValueError("Unexpected number of columns in given BED file")
    region_names = bed_df_copy[3]
    duplicated_region_names = region_names[region_names.duplicated()]
    if duplicated_region_names.any():
        logging.error(
            f"Generated duplicate region names from BED file: {list(duplicated_region_names)}. This could indicate that there are multiple regions in the bed file with the same coordinates."
        )
        sys.exit(1)
    return bed_df_copy


class Region:
    def __init__(
        self,
        genome_depth_per_base: pd.DataFrame,
        start: int,
        end: int,
        name: str,
        sample: str,
        coverage_thresholds: list[float] = [],
    ):
        self.region_depth = self._get_region_depth_per_base(
            genome_depth_per_base, start, end
        )
        self.start = start
        self.end = end
        self.name = name
        self.sample = sample
        self.summary_stats = self._get_summary_stats()
        self.coverage = self.get_percent_coverage(coverage_thresholds)

    def _get_summary_stats(self) -> dict:
        depth_per_base = self.region_depth
        return {
            "length": len(depth_per_base["Cov"]),
            "missing_sites": sum(depth_per_base["Cov"] == 0),
            "depth_median": np.median(depth_per_base["Cov"]),
            "depth_mean": np.mean(depth_per_base["Cov"]),
            "depth_min": np.min(depth_per_base["Cov"]),
            "depth_max": np.max(depth_per_base["Cov"]),
        }

    def get_percent_coverage(self, thresholds: list = []) -> dict:
        if len(thresholds) == 0:
            return {}
        coverage_at_thresholds = {}
        for threshold in thresholds:
            coverage_at_thresholds[f"cov_perc_{threshold}x"] = (
                self._get_percent_coverage_at_threshold(threshold)
            )
        return coverage_at_thresholds

    def _get_percent_coverage_at_threshold(self, threshold: float) -> float:
        """
        Compute percentage of bases covered at given (coverage) threshold
        """
        return (
            sum(self.region_depth["Cov"] >= threshold) / len(self.region_depth)
        ) * 100

    @staticmethod
    def _get_region_depth_per_base(
        genome_depth_per_base: pd.DataFrame, start: int, end: int
    ) -> pd.DataFrame:
        return genome_depth_per_base.loc[
            (genome_depth_per_base["Pos"] >= int(start) - 1)
            & (genome_depth_per_base["Pos"] <= int(end) - 1)
        ]

    def to_dict(self):
        return {
            k: v
            for (k, v) in vars(self).items()
            if not k.startswith("_") and not callable(v)
        }

    def to_row(self):
        row = self.to_dict()
        keys_to_remove = ["region_depth", "summary_stats", "coverage"]
        for k in keys_to_remove:
            del row[k]
        row.update(self.summary_stats)
        row.update(self.coverage)
        return row


## OUTPUT
def generate_coverage_summary_rows(
    depth_df: pd.DataFrame, bed_df: pd.DataFrame, coverage_thresholds: list[float]
) -> list[dict]:
    rows = []
    for row_index in range(0, len(bed_df.index)):
        region = Region(
            depth_df,
            bed_df.iloc[row_index, 1],
            bed_df.iloc[row_index, 2],
            bed_df.iloc[row_index, 3],
            sample_name,
            coverage_thresholds,
        )
        rows.append(region.to_row())
    return rows


def dataframe_from_rows(
    rows: list[dict], header_override: dict = None, header_order: list = []
) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if header_override is not None:
        df = df.rename(columns=header_override)
    if len(header_order) == 0:
        header_order = sorted(df.columns)
    return df[header_order]


def make_output_dir(output: Path) -> None:
    if not output.is_dir() and output.exists():
        logging.error(
            f"Given output {output} exists but is not a directory. Please supply a valid directory path."
        )
        sys.exit(1)
    else:
        output.mkdir(parents=True, exist_ok=True)


def write_output_csv(
    rows: list[dict],
    thresholds: list[str],
    output_dir: Path,
    output_filename: str = "coverage_summary.tsv",
    precision: int = 1,
) -> None:
    if len(rows) == 0:
        raise ValueError(
            "Unexpected number of rows generated from given BED file and samtools depth output."
        )
    header = [
        "sample",
        "name",
        "start",
        "end",
        "length",
        "missing_sites",
        "depth_median",
        "depth_mean",
        "depth_min",
        "depth_max",
    ]
    coverage_threshold_headers = [f"cov_perc_{threshold}x" for threshold in thresholds]
    header.extend(coverage_threshold_headers)
    dfout = dataframe_from_rows(rows, header_order=header)
    dfout.to_csv(
        output_dir / output_filename,
        sep="\t",
        index=False,
        header=True,
        float_format=f"%.{precision}f",
    )


if __name__ == "__main__":
    args = parse_args()

    # Make output dir
    make_output_dir(args.output)

    # Get sample name
    if not args.sample_name:
        sample_name = str(args.samtools_depth.stem).removesuffix("_samtools_depth")
    else:
        sample_name = args.sample_name

    # Parse depth file
    depth_df = parse_depth_data(args.samtools_depth)

    # Parse bed file
    bed_df = parse_bed_file(args.bed)

    # Analyse regions defined in bed file
    rows = generate_coverage_summary_rows(depth_df, bed_df, args.coverage_thresholds)

    # Output to tsv
    write_output_csv(
        rows,
        args.coverage_thresholds,
        output_dir=args.output,
        output_filename=f"{sample_name}_coverage_summary.tsv",
    )
