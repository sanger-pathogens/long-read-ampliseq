#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd
# import numpy as np
import plotly.graph_objs as go
# from plotly.subplots import make_subplots


def parse_args():
    parser = argparse.ArgumentParser(description="A script to generate coverage plots from tabulated data")

    parser.add_argument("--coverage_summary_dir", "-d", help="Directory containing coverage summary TSV files")
    return parser.parse_args()

def parse_tsvs(tsv_files: list[Path]) -> pd.DataFrame:
    return pd.concat((pd.read_csv(f, sep='\t') for f in tsv_files), ignore_index=True)

def extract_matrix(df: pd.DataFrame, cov_col="cov_perc_25.0x") -> pd.DataFrame:
    return pd.pivot(df, index="name", columns="sample", values=cov_col)

def get_cov_columns(df: pd.DataFrame) -> list:
    return [col for col in df.columns if col.startswith("cov_perc")] 

def plot_heatmap(df: pd.DataFrame, cov="25.0x") -> None:
    fig = go.Figure(
        data=go.Heatmap(
            z=df.to_numpy(),
            x=df.columns,
            y=df.index,
            colorscale='viridis'
        )
    )

    fig.update_layout(
        title=f'Coverage at {cov}'
    )
    fig.update_xaxes(title_text="Sample")
    fig.update_yaxes(title_text="Locus")
    
    fig.show()

def main():
    args = parse_args()

    # Get all TSV files in the directory
    tsv_files = list(Path(args.coverage_summary_dir).glob("*.tsv"))

    # Exit if no tsvs
    if len(tsv_files) == 0:
        sys.exit(1, "No TSV files found in the directory.")

    df = parse_tsvs(tsv_files)
    cov_cols = get_cov_columns(df)
    for col in cov_cols:
        cov_str = col.split("_")[-1]
        sliced_df = extract_matrix(df, cov_col=col)
        plot_heatmap(sliced_df, cov=cov_str)

if __name__ == "__main__":
    main()