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
    parser.add_argument("-s", type=str, dest="samtools_depth",help="Specify path to output from `samtools depth`")
    parser.add_argument("-b", type=Path, dest="bam_file",help="Specify input bam file")
    parser.add_argument("-w", dest="wgs", action="store_true", help="Calculate stats for Whole Genome" )
    parser.add_argument("-r", type=str, dest="region", help="Comma separated region positions (<start>,<end>,<name>)" )
    parser.add_argument("-c", type=str, dest="bed", help="BED file (TSV) defining regions (<name>\t<start>\t<end>)" )
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


def get_covwindow(depth_per_base, start, end, region_name):
    """
    Extract depth for a given region and calculate summary coverage stats
    """
    region_depth = get_region_depth_per_base(depth_per_base, start, end)
    region_length = len(region_depth['Cov'])
    missing_sites = sum(region_depth['Cov'] == 0)
    region_median = np.median(region_depth['Cov'])
    region_mean = np.mean(region_depth['Cov'])
    mincov = np.min(region_depth['Cov'])
    maxcov = np.max(region_depth['Cov'])
    perc1x = round((len(list(region_depth['Cov'][region_depth['Cov']>=1]))/region_length)*100,1)
    perc5x = round((len(list(region_depth['Cov'][region_depth['Cov']>=5]))/region_length)*100,1)
    perc8x = round((len(list(region_depth['Cov'][region_depth['Cov']>=8]))/region_length)*100,1)
    perc20x = round((len(list(region_depth['Cov'][region_depth['Cov']>=20]))/region_length)*100,1)
    perc100x = round((len(list(region_depth['Cov'][region_depth['Cov']>=100]))/region_length)*100,1)
    perc250x = round((len(list(region_depth['Cov'][region_depth['Cov']>=250]))/region_length)*100,1)

    return [mysamplename, region_name, str(start), str(end), str(region_length), str(missing_sites), str(region_median), str(region_mean), str(mincov), str(maxcov), str(perc1x), str(perc5x), str(perc8x), str(perc20x), str(perc100x), str(perc250x)]


def get_percent_coverage_at_threshold(depth_per_base: pd.Series, threshold: float) -> float:
    """
    Compute percentage of bases covered at given (coverage) threshold and round to 1 decimal place
    """
    percent_coverage = sum(depth_per_base >= threshold) / len(depth_per_base) * 100
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
        coverage_at_thresholds[f"perc{threshold}x"] = get_percent_coverage_at_threshold(depth_per_base, threshold)
    return coverage_at_thresholds


if __name__ == "__main__":
    args = parse_args()
    bam_file = args.bam_file

    # Get sequence name 
    mysamplename = bam_file.stem

    mydepth = pd.read_csv(args.samtools_depth, sep='\t',index_col=0, header=None,names=['Ref','Pos','Cov'])

    # Anything that isn't covered (is na) make 0
    mydepth[pd.isnull(mydepth['Cov'])] = 0

    # Calculate genomewide median and mean coverage (median is probably more reliable here)
    full_median = np.median(mydepth['Cov'])
    full_mean = np.mean(mydepth['Cov'])
    wgs_length = len(list(mydepth['Pos']))



    header = ["Sample", "Region", "StartPos", "EndPos", "Length", "Missing_Sites", "Median_Cov", "Mean_Cov", "Min_Cov", "Max_Cov", "cov1x.perc", "cov5x.perc", "cov8x.perc", "cov20x.perc", "cov100x.perc", "cov250x.perc"]
    dfout = pd.DataFrame(columns = header)

    # Run function on a whole genome
    if args.wgs :
        region_cov_summary = get_covwindow(mydepth, "1", wgs_length, "WGS")
        dfout = pd.concat([pd.DataFrame([region_cov_summary], columns=header)], ignore_index=True)


    # Run function on a single region
    if args.region :
        myregion = [s.strip() for s in args.region.split(",")]
        region_cov_summary = get_covwindow(mydepth, myregion[0], myregion[1], myregion[2])
        dfout = pd.concat([pd.DataFrame([region_cov_summary], columns=header)], ignore_index=True)



    # Run function on a csv list of regions
    if args.bed : 
        bed_file = pd.read_csv(args.bed, header=None, sep='\t')
        print(bed_file)
        for myline in range(0, len(bed_file.index)):
            region_cov_summary = get_covwindow(mydepth, bed_file.iloc[myline,1], bed_file.iloc[myline,2], bed_file.iloc[myline,0])
            dfout = pd.concat([dfout, pd.DataFrame([region_cov_summary], columns=header)], ignore_index=True)


    #print(dfout)
    dfout.to_csv(mysamplename + ".depth.tsv", sep="\t", index=False, header=True)

