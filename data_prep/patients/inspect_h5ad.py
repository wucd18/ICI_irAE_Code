#!/usr/bin/env python3
"""Audit H5AD structure and export cell metadata without loading expression into RAM."""

from __future__ import annotations

import argparse
import json
import pathlib

import anndata as ad
import pandas as pd


def safe_values(series: pd.Series, limit: int = 30) -> list[str] | None:
    count = int(series.nunique(dropna=False))
    if count > limit:
        return None
    values = series.astype("string").fillna("<NA>").drop_duplicates().tolist()
    return [str(value) for value in values]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("summary_json")
    parser.add_argument("obs_parquet")
    args = parser.parse_args()

    input_path = pathlib.Path(args.input).resolve()
    dataset = ad.read_h5ad(input_path, backed="r")
    obs = dataset.obs.copy()
    var = dataset.var.copy()

    columns = []
    for column in obs.columns:
        series = obs[column]
        columns.append(
            {
                "column": str(column),
                "dtype": str(series.dtype),
                "n_unique": int(series.nunique(dropna=False)),
                "n_missing": int(series.isna().sum()),
                "values_if_low_cardinality": safe_values(series),
            }
        )

    summary = {
        "input": str(input_path),
        "n_obs": int(dataset.n_obs),
        "n_vars": int(dataset.n_vars),
        "obs_columns": columns,
        "var_columns": [str(column) for column in var.columns],
        "layers": list(dataset.layers.keys()),
        "obsm": list(dataset.obsm.keys()),
        "obsp": list(dataset.obsp.keys()),
        "varm": list(dataset.varm.keys()),
        "uns": list(dataset.uns.keys()),
        "raw_present": dataset.raw is not None,
        "x_type": type(dataset.X).__name__,
    }

    pathlib.Path(args.summary_json).write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    obs.reset_index(names="cell_barcode").to_parquet(args.obs_parquet, index=False)
    dataset.file.close()
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
