#!/usr/bin/env python
import sys
import io
import subprocess
import os
import pandas as pd
import numpy as np
import re
import argparse
import csv
import logging

from pathlib import Path

# Copyright Mat Beale, Wellcome Sanger Institute, July 2022


def parse_args():
    parser = argparse.ArgumentParser(
        description="Takes a bam file and a list of regions, and uses samtools depth to determine the mean depth per region, the normalised depth (useful for CNV detection), and the number of missing (zero coverage) sites per gene. Note that this script currently assumes a single reference chromosome."
    )
    parser.add_argument("-s", type=Path, dest="samtools_depth",help="Specify path to output from `samtools depth`")
    parser.add_argument("-b", type=Path, dest="bam_file",help="Specify input bam file")
    parser.add_argument("-c", type=Path, dest="bed", help="BED file (TSV) defining regions (<chrom>\t<start>\t<end>\t<name>)" )
    parser.add_argument("-t", type=parse_coverage_threshold, dest="coverage_threshold", help="Comma separated list of coverage thresholds")
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
            "median": np.median(depth_per_base['Cov']),
            "mean": np.mean(depth_per_base['Cov']),
            "min": np.min(depth_per_base['Cov']),
            "max": np.max(depth_per_base['Cov']),
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
    
    



def get_covwindow(depth_per_base, start, end, region_name, sample, thresholds):
    """
    Extract depth for a given region and calculate summary coverage stats
    """
    region_depth = get_region_depth_per_base(depth_per_base, start, end)
    summary_stats = compute_region_summary_stats(region_depth)
    coverage_at_thresholds = compute_region_percent_coverage_over_thresholds(depth_per_base, thresholds)

    summary_stat_order = ["length", "missing_sites", "median", "mean", "min", "max"]
    coverage_order = sorted(coverage_at_thresholds.keys())

    output_region_stats = [sample, region_name, start, end]
    output_region_stats.extend([summary_stats[k] for k in summary_stat_order])
    output_region_stats.extend([coverage_at_thresholds[k] for k in coverage_order])

    return output_region_stats


def get_percent_coverage_at_threshold(depth_per_base: pd.Series, threshold: float) -> float:
    """
    Compute percentage of bases covered at given (coverage) threshold and round to 1 decimal place
    """
    print(depth_per_base['Cov'] >= threshold)
    percent_coverage = (sum(depth_per_base['Cov'] >= threshold) / len(depth_per_base)) * 100
    print(percent_coverage)
    exit
    return round(percent_coverage, 1)


def get_region_depth_per_base(genome_depth_per_base: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    return genome_depth_per_base.loc[(genome_depth_per_base['Pos'] >= int(start)-1) & (genome_depth_per_base['Pos'] <= int(end)-1)]


def compute_region_summary_stats(depth_per_base: pd.DataFrame) -> dict:
    return {
        "length": len(depth_per_base['Cov']),
        "missing_sites": sum(depth_per_base['Cov'] == 0),
        "median": np.median(depth_per_base['Cov']),
        "mean": np.mean(depth_per_base['Cov']),
        "min": np.min(depth_per_base['Cov']),
        "max": np.max(depth_per_base['Cov']),
    }


def compute_region_percent_coverage_over_thresholds(depth_per_base: pd.DataFrame, thresholds: list) -> dict:
    coverage_at_thresholds = {}
    for threshold in thresholds:
        coverage_at_thresholds[f"cov_perc_{threshold}x"] = get_percent_coverage_at_threshold(depth_per_base, threshold)
    return coverage_at_thresholds


def dataframe_from_rows(rows: list[dict], header_override: dict = None, header_order: list = []) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if header_override is not None:
        df = df.rename(columns=header_override)
    if len(header_order) == 0:
        header_order = sorted(df.columns)
    return df[header_order]
    

if __name__ == "__main__":
    args = parse_args()
    bam_file = args.bam_file

    # Get sequence name 
    sample_name = bam_file.stem

    mydepth = pd.read_csv(args.samtools_depth, sep='\t', index_col=0, header=None, names=['Ref','Pos','Cov'])

    # Anything that isn't covered (is na) make 0
    mydepth[pd.isnull(mydepth['Cov'])] = 0

    # Run function on a csv list of regions
    header_override = {"median": "depth_median", "mean": "depth_mean", "min": "depth_min", "max": "depth_max"}
    header = ["sample", "name", "start", "end", "length", "missing_sites", "depth_median", "depth_mean", "depth_min", "depth_max"]
    coverage_threshold_headers = [f"cov_perc_{threshold}x" for threshold in args.coverage_threshold]
    header.extend(coverage_threshold_headers)
    if args.bed : 
        bed_file = pd.read_csv(args.bed, header=None, sep='\t')
        rows = []
        for myline in range(0, len(bed_file.index)):
            region = Region(mydepth, bed_file.iloc[myline,1], bed_file.iloc[myline,2], bed_file.iloc[myline,3], sample_name, args.coverage_threshold)
            rows.append(region.to_row())
        dfout = dataframe_from_rows(rows, header_override, header_order=header)
            # region_cov_summary = get_covwindow(mydepth, bed_file.iloc[myline,1], bed_file.iloc[myline,2], bed_file.iloc[myline,0], sample_name, args.coverage_threshold)
            # print(region_cov_summary)
            # dfout = pd.concat([dfout, pd.DataFrame([region_cov_summary], columns=header)], ignore_index=True)

    dfout.to_csv(sample_name + ".depth.tsv", sep="\t", index=False, header=True)

