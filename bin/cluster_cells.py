#!/usr/bin/env python3
import argparse
import anndata as ad
import pandas as pd
import scanpy as sc


def _sanitize_strings_for_compat(adata):
    """Write standard string-array datasets so R/zellkonverter (<0.11) can read them."""
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
    parser = argparse.ArgumentParser(description="Cluster cells to build a global expression baseline.")
    parser.add_argument("--h5ad", required=True, help="Input h5ad file")
    parser.add_argument("--resolution", type=float, default=1.0, help="Leiden clustering resolution")
    parser.add_argument("--out_annotations", required=True, help="Output TSV for InferCNV")
    parser.add_argument("--out_h5ad", required=True, help="Output annotated h5ad file")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)

    # Honor upstream QC if present (e.g. scrnaseq SCANPY_QC writes passing_qc
    # but does not subset).
    if "passing_qc" in adata.obs.columns:
        adata = adata[adata.obs["passing_qc"]].copy()

    bdata = adata.copy()

    sc.pp.normalize_total(bdata, target_sum=1e4)
    sc.pp.log1p(bdata)
    sc.pp.highly_variable_genes(bdata, n_top_genes=2000)
    sc.pp.pca(bdata)
    sc.pp.neighbors(bdata)
    sc.tl.leiden(bdata, resolution=args.resolution)

    adata.obs['leiden'] = 'cluster_' + bdata.obs['leiden'].astype(str)
    adata.obs[['leiden']].to_csv(args.out_annotations, sep="\t", header=False)
    _sanitize_strings_for_compat(adata)
    adata.write_h5ad(args.out_h5ad)

if __name__ == "__main__":
    main()
