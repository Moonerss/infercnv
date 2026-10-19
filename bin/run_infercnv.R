#!/usr/bin/env Rscript
library(optparse)
library(infercnv)
library(zellkonverter)
library(SummarizedExperiment)
library(Matrix)

option_list <- list(
    make_option(c("--h5ad"), type="character", default=NULL, help="Input H5AD file"),
    make_option(c("--annotations"), type="character", default=NULL, help="Barcode to cluster annotations"),
    make_option(c("--gene_order"), type="character", default=NULL, help="Genomic position file"),
    make_option(c("--cutoff"), type="double", default=0.1, help="Count cutoff threshold"),
    make_option(c("--threads"), type="integer", default=4L, help="Number of threads"),
    make_option(c("--out_dir"), type="character", default=".", help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$h5ad) || is.null(opt$annotations) || is.null(opt$gene_order)) {
    print_help(opt_parser)
    stop("Missing required arguments.", call.=FALSE)
}

# anndata >= 0.11 + pandas >= 3.0 writes string indices/columns with
# encoding-type='nullable-string-array', which the older anndata pinned in
# zellkonverter's basilisk environment (< 0.11) cannot read. Downgrade any
# nullable-string-array groups to their underlying string-array 'values'
# dataset in a temporary copy before calling readH5AD().
compat_h5ad <- basilisk::basiliskRun(
    env = zellkonverter::zellkonverterAnnDataEnv(),
    fun = function(in_path) {
        reticulate::py_run_string("
import shutil
import h5py

def _has_nullable_strings(grp):
    for k in grp.keys():
        obj = grp[k]
        if isinstance(obj, h5py.Group):
            if obj.attrs.get('encoding-type') == 'nullable-string-array':
                return True
            if _has_nullable_strings(obj):
                return True
    return False

def _downgrade_nullable_strings(grp):
    for k in list(grp.keys()):
        obj = grp[k]
        if isinstance(obj, h5py.Group):
            if obj.attrs.get('encoding-type') == 'nullable-string-array' and 'values' in obj:
                tmp_key = k + '__compat_values__'
                grp.move(k + '/values', tmp_key)
                del grp[k]
                grp.move(tmp_key, k)
            else:
                _downgrade_nullable_strings(obj)

def ensure_compat_h5ad(in_path, out_path):
    with h5py.File(in_path, 'r') as f:
        if not _has_nullable_strings(f):
            return in_path
    shutil.copyfile(in_path, out_path)
    with h5py.File(out_path, 'r+') as f:
        _downgrade_nullable_strings(f)
    return out_path
")
        reticulate::py$ensure_compat_h5ad(in_path, "compat_input.h5ad")
    },
    in_path = opt$h5ad
)

adata <- tryCatch(
    readH5AD(compat_h5ad),
    finally = {
        if (compat_h5ad != opt$h5ad && file.exists(compat_h5ad)) {
            unlink(compat_h5ad)
        }
    }
)

# zellkonverter names the primary assay based on adata.uns['X_name'] (often
# "X" or "counts") — grab the first assay rather than a hardcoded name.
expression <- as.matrix(assays(adata)[[1]])

infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = expression,
    annotations_file  = opt$annotations,
    gene_order_file   = opt$gene_order,
    ref_group_names   = NULL
)

infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff            = opt$cutoff,
    out_dir           = opt$out_dir,
    cluster_by_groups = TRUE,
    denoise           = TRUE,
    HMM               = TRUE,
    analysis_mode     = "subclusters",
    no_plot           = FALSE,
    num_threads       = opt$threads
)
