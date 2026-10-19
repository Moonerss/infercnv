#!/usr/bin/env python3
import argparse
import anndata as ad
import pandas as pd
import scanpy as sc


def _sanitize_strings_for_compat(adata):
    """Write standard string-array datasets for downstream compatibility."""
    if hasattr(ad, "settings") and hasattr(ad.settings, "allow_write_nullable_strings"):
        ad.settings.allow_write_nullable_strings = False
    frames = [adata.obs, adata.var]
    if adata.raw is not None:
        frames.append(adata.raw.var)
    for df in frames:
        df.index = df.index.astype(object)
        for col in df.columns:
            if isinstance(df[col].dtype, pd.StringDtype):
                df[col] = df[col].fillna("").astype(object)
            elif isinstance(df[col].dtype, pd.CategoricalDtype) and isinstance(
                df[col].cat.categories.dtype, pd.StringDtype
            ):
                df[col] = df[col].cat.rename_categories(
                    df[col].cat.categories.astype(object)
                )


def main():
    parser = argparse.ArgumentParser(description="Concatenate an array of h5ad files.")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--out_h5ad", required=True)
    args = parser.parse_args()

    adatas = [sc.read_h5ad(f) for f in args.inputs]

    # index_unique avoids barcode collisions across samples (10x barcodes can
    # repeat across libraries).
    combined = ad.concat(adatas, merge="same", index_unique="-")
    _sanitize_strings_for_compat(combined)
    combined.write_h5ad(args.out_h5ad)


if __name__ == "__main__":
    main()
