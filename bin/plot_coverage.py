#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objs as go

DEFAULT_HEATMAP_FILENAME = "coverage_heatmap.html"

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def parse_args():
    parser = argparse.ArgumentParser(description="A script to generate coverage plots from tabulated data")

    parser.add_argument("coverage_summary_dir", help="Directory containing coverage summary TSV files", type=Path)
    parser.add_argument("--interactive", "-i", help="Interactive mode to view plots, instead of writing to a file", action="store_true")
    parser.add_argument("--output", "-o", help="Path to output directory", type=Path, default=".")
    return parser.parse_args()

def validate_coverage_summary_dir(coverage_summary_dir: Path) -> Path:
    if not coverage_summary_dir.is_dir():
        logging.error(f"Given coverage summary directory path '{coverage_summary_dir}' does not exist or is not a directory.")
        sys.exit(1)
    return coverage_summary_dir

def validate_output_arg(output: Path) -> Path:
    if output.exists() and not output.is_dir():
        logging.error(f"Path provided to --output exists, but is not a directory: '{output}'")
        sys.exit(1)
    return output

def extract_tsvs(directory: Path) -> list[Path]:
    """Get paths to all TSV files in the given directory"""
    tsv_files = list(directory.glob("*.tsv"))
    if len(tsv_files) == 0:
        logging.error(f"No TSV files found in directory: '{directory}'")
        sys.exit(1)
    return tsv_files

def parse_tsvs(tsv_files: list[Path]) -> pd.DataFrame:
    return pd.concat((pd.read_csv(f, sep='\t') for f in tsv_files), ignore_index=True)

def extract_matrix(df: pd.DataFrame, cov_col="cov_perc_25.0x") -> pd.DataFrame:
    return pd.pivot(df, index="name", columns="sample", values=cov_col)

def get_cov_columns(df: pd.DataFrame) -> list:
    return [col for col in df.columns if col.startswith("cov_perc")]

def get_cov_level(cov_col: str) -> float:
    return float(cov_col.split("_")[-1].rstrip("x"))

def get_cov_levels(cov_cols: list) -> list:
    return [get_cov_level(col) for col in cov_cols]


def add_heatmap_traces(fig: go.Figure, df: pd.DataFrame, cov_cols: list) -> go.Figure:
    cov_levels = get_cov_levels(cov_cols)
    for cov, col in zip(cov_levels, cov_cols):
        sliced_df = extract_matrix(df, cov_col=col)
        fig.add_trace(
            go.Heatmap(
                z=sliced_df.to_numpy(),
                x=sliced_df.columns,
                y=sliced_df.index,
                colorscale='viridis',
                name=f"cov = {cov}x",
                zmin=0,
                zmax=100
            )
        )
    return fig

def create_cov_slider(cov_levels: list[float]) -> dict:
    steps = []
    for i, cov in enumerate(cov_levels):
        step = dict(
            method="restyle",
            args=[
                {"visible": [False] * len(cov_levels), "title": "Coverage: " + str(cov)}
            ],
            label=str(cov)
        )
        step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
        steps.append(step)
    
    return dict(
        active=0,
        currentvalue={"prefix": "Coverage: "},
        pad={"t": 100},
        steps=steps
    )

def format_heatmap_figure(fig: go.Figure, sliders: list) -> go.Figure:
    # Need to ensure that only first trace is visible initially
    for i, trace in enumerate(fig.data):
        if i == 0:
            trace.visible = True
        else:
            trace.visible = False

    fig.update_layout(
        title=f'Coverage',
        sliders=sliders
    )
    fig.update_traces(colorbar_title_text="% bases")
    fig.update_xaxes(title_text="Sample")
    fig.update_yaxes(title_text="Locus")
    
    return fig

def main():
    args = parse_args()

    # Validate args
    coverage_summary_dir = validate_coverage_summary_dir(args.coverage_summary_dir)
    output = validate_output_arg(args.output)
    
    # Parse
    tsv_files = extract_tsvs(coverage_summary_dir)

    df = parse_tsvs(tsv_files)
    cov_cols = get_cov_columns(df)
    cov_levels = get_cov_levels(cov_cols)
    
    # Create Heatmap
    fig = go.Figure()
    cov_slider = create_cov_slider(cov_levels)
    fig = add_heatmap_traces(fig, df, cov_cols)
    fig = format_heatmap_figure(fig, sliders=[cov_slider])

    # Output
    if args.interactive:
        fig.show()
        if output:
            logging.warn(f"Running in interactive mode. An output argument was provided, but was ignored.")
    else:
        output.mkdir(parents=True, exist_ok=True)
        fig.write_html(output / DEFAULT_HEATMAP_FILENAME)


if __name__ == "__main__":
    main()