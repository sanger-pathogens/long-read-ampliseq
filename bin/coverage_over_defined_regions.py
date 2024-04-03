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
    parser.add_argument("--samtools-depth", "-s", type=parse_samtools_depth, required=True, dest="samtools_depth",help="Path to `samtools depth` output file")
    parser.add_argument("--bed", "-b", type=parse_bed_file, required=True, dest="bed", help=r"Path to BED file (TSV) defining regions (<chrom>\t<start>\t<end>\t<name>)" )
    parser.add_argument("--threshold", "-t", type=parse_coverage_threshold, required=True, dest="coverage_threshold", help="Comma separated list of coverage thresholds")
    parser.add_argument("--sample-name", "-n", type=str, required=True, dest="sample_name", help="Sample name")
    parser.set_defaults(feature=False)

    return parser.parse_args()

def parse_coverage_threshold(coverage_threshold_arg: str) -> list[float]:
    valid_thresholds = []
    coverage_thresholds = coverage_threshold_arg.split(",")
    for threshold in coverage_thresholds:
        try:
            valid_threshold = float(threshold)
        except:
            logging.error(f"Given threshold '{threshold}' is not a valid threshold value.")
            sys.exit(1)
        valid_thresholds.append(valid_threshold)
    return valid_thresholds

def parse_bed_file(path: str) -> pd.DataFrame:
    bed_path = Path(path)
    if not bed_path.is_file():
        logging.error(f"Given path '{path}' is not a file.")
        sys.exit(1)
    else:
        df = pd.read_csv(path, header=None, sep='\t')
    if len(df.columns) != 4:
        logging.error(f"Given BED file '{bed_path}' is not a valid BED4 file (please ensure the file contains 4 columns).")
        sys.exit(1)
    return df

def parse_samtools_depth(path: str) -> Path:
    samtools_depth_path = Path(path)
    if not samtools_depth_path.is_file():
        logging.error(f"Given samtools depth output file '{path}' is not a valid file.")
        sys.exit(1)
    return samtools_depth_path


class Region:
    def __init__(self, genome_depth_per_base: pd.DataFrame, start: int, end: int, name: str, sample: str, coverage_thresholds:list[float] = []):
        self.region_depth = self._get_region_depth_per_base(genome_depth_per_base, start, end)
        self.start = start
        self.end = end
        self.name = name
        self.sample = sample
        self.summary_stats = self._get_summary_stats()
        self.coverage = self.get_percent_coverage(coverage_thresholds)

    def _get_summary_stats(self) -> dict:
        depth_per_base = self.region_depth
        return {
            "length": len(depth_per_base['Cov']),
            "missing_sites": sum(depth_per_base['Cov'] == 0),
            "depth_median": np.median(depth_per_base['Cov']),
            "depth_mean": np.mean(depth_per_base['Cov']),
            "depth_min": np.min(depth_per_base['Cov']),
            "depth_max": np.max(depth_per_base['Cov']),
        }
    
    def get_percent_coverage(self, thresholds: list = []) -> dict:
        if len(thresholds) == 0:
            return {}
        coverage_at_thresholds = {}
        for threshold in thresholds:
            coverage_at_thresholds[f"cov_perc_{threshold}x"] = self._get_percent_coverage_at_threshold(threshold)
        return coverage_at_thresholds

    def _get_percent_coverage_at_threshold(self, threshold: float) -> float:
        """
        Compute percentage of bases covered at given (coverage) threshold
        """
        return (sum(self.region_depth['Cov'] >= threshold) / len(self.region_depth)) * 100
    
    @staticmethod
    def _get_region_depth_per_base(genome_depth_per_base: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
        return genome_depth_per_base.loc[(genome_depth_per_base['Pos'] >= int(start)-1) & (genome_depth_per_base['Pos'] <= int(end)-1)]

    def to_dict(self):
        return {k: v for (k, v) in vars(self).items() if not k.startswith("_") and not callable(v)}

    def to_row(self):
        row = self.to_dict()
        keys_to_remove = ["region_depth", "summary_stats", "coverage"]
        for k in keys_to_remove:
            del row[k]
        row.update(self.summary_stats)
        row.update(self.coverage)
        return row


def dataframe_from_rows(rows: list[dict], header_override: dict = None, header_order: list = []) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if header_override is not None:
        df = df.rename(columns=header_override)
    if len(header_order) == 0:
        header_order = sorted(df.columns)
    return df[header_order]
    

if __name__ == "__main__":
    args = parse_args()

    # Get sample name
    if not args.sample_name:
        sample_name = str(args.samtools_depth.stem).removesuffix("_samtools_depth")
    else:
        sample_name = args.sample_name

    # Parse depth file
    mydepth = pd.read_csv(args.samtools_depth, sep='\t', index_col=0, header=None, names=['Ref','Pos','Cov'])
    # Anything that isn't covered (is na) make 0
    mydepth[pd.isnull(mydepth['Cov'])] = 0

    # Analyse regions defined in bed file
    bed_file = args.bed
    rows = []
    for row_index in range(0, len(bed_file.index)):
        region = Region(mydepth, bed_file.iloc[row_index, 1], bed_file.iloc[row_index, 2], bed_file.iloc[row_index, 3], sample_name, args.coverage_threshold)
        rows.append(region.to_row())

    # Output to tsv
    if rows:
        header = ["sample", "name", "start", "end", "length", "missing_sites", "depth_median", "depth_mean", "depth_min", "depth_max"]
        coverage_threshold_headers = [f"cov_perc_{threshold}x" for threshold in args.coverage_threshold]
        header.extend(coverage_threshold_headers)
        dfout = dataframe_from_rows(rows, header_order=header)
        dfout.to_csv(sample_name + ".depth.tsv", sep="\t", index=False, header=True, float_format="%.1f")

