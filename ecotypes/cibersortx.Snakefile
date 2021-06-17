# ==================================
# Dan Roden
# CREATED: 13/07/2020
# ----------------------------------
# Snakemake file for running CIBERSORTx decon
# NOTE: actual CIBERSORTx runs are done via a separate shell script: "run_docker_cibersortx.sh"
# ------------------------------------------------------------------------------
# ==========
# CONFIG
# ==========
# No config.json used
#
# ==================================
from snakemake.utils import R
from pprint import pprint
import itertools

SUBSAMPLE_FRACTION=[0.15]

# Sample ID Mappings
ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020={"Jul2020": "analysis/Jul2020_updated_objectfile/Rdata_final/01_miniatlas_final_object.Rdata"}

NORMAL_EPITHELIAL_CELL_CALLS=["with_normal_epithelial_cells"]
CELL_TYPES=["celltype_subset"]
CANCER_CELL_ANNOTATIONS=["SCTyper_Dec"]
CYCLING_CELL_ANNOTATIONS=["combined"]

# ==================================
rule all:
    input:
        # generate_cell_barcode_lists_from_seurat_metadata
        expand("analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/barcode_list.txt", analysis=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020.keys(), normal_epithelial_cells=NORMAL_EPITHELIAL_CELL_CALLS, cell_type=CELL_TYPES, cancer_cell_annotation=CANCER_CELL_ANNOTATIONS, cycling_cell_annotation=CYCLING_CELL_ANNOTATIONS),
        expand("analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/{files}.txt", analysis=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020.keys(), normal_epithelial_cells=NORMAL_EPITHELIAL_CELL_CALLS, cell_type=CELL_TYPES, cancer_cell_annotation=CANCER_CELL_ANNOTATIONS, cycling_cell_annotation=CYCLING_CELL_ANNOTATIONS, subsample_fraction=SUBSAMPLE_FRACTION, files=["counts_sparse_matrix", "cells", "genes"]),
        # convert_sparse_matrix_to_dense
        expand("analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/counts_dense_matrix.txt", analysis=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020.keys(), normal_epithelial_cells=NORMAL_EPITHELIAL_CELL_CALLS, cell_type=CELL_TYPES, cancer_cell_annotation=CANCER_CELL_ANNOTATIONS, cycling_cell_annotation=CYCLING_CELL_ANNOTATIONS, subsample_fraction=SUBSAMPLE_FRACTION),
        # generate_pseudobulk_mixture_file
        expand("analysis/CIBERSORTx/{analysis}/pseudobulk_mixture_matrix/sum_matching_pseudobulk_mixture_matrix.txt", analysis=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020.keys()),
        expand("analysis/CIBERSORTx/{analysis}/pseudobulk_mixture_matrix/avg_matching_pseudobulk_mixture_matrix.txt", analysis=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020.keys())


# ==================================
# Generate cell barcode lists (from seurat meta-data)
# ----------------------------------
# Extract barcodes of cell subsets based on a meta-data column
# ==================================
rule generate_cell_barcode_lists_July2020:
    input:
        seurat=lambda wildcards: expand("{analysis_data}", analysis_data=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020[wildcards.analysis])
    output:
        barcode_list="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/barcode_list.txt"
    params:
        cores="1",
        hvmem="h_vmem=100G",
        memrequested="mem_requested=100G",
        tmprequested="tmp_requested=10G",
        tmpfree="tmpfree=10G",
        jobname="generate_cell_barcode_lists_Sept2019"
    script:
        "code/celltype_deconvolution/generate_cell_barcode_lists.snakemake.R"

# ==================================
# extract_sparse_matrix_for_signature_matrix
# ----------------------------------
# Extract sparse count matrix for signature matrix generation
# ==================================
rule extract_sparse_matrix_for_celltype_signatures:
    input:
        seurat=lambda wildcards: expand("{analysis_data}", analysis_data=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020[wildcards.analysis]),
        barcode_list="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/barcode_list.txt"
    output:
        sparse_matrix="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/counts_sparse_matrix.txt",
        cells="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/cells.txt",
        genes="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/genes.txt"
    params:
        cores="1",
        hvmem="h_vmem=200G",
        memrequested="mem_requested=200G",
        tmprequested="tmp_requested=10G",
        tmpfree="tmpfree=10G",
        jobname="extract_sparse_matrix_for_celltype_signatures"
    script:
        "code/celltype_deconvolution/extract_sparse_matrix_for_celltype_signatures_from_barcode_list.snakemake.R"

# ==================================
# convert_sparse_matrix_to_dense
# ----------------------------------
# Convert sparse count matrix to dense form for signature matrix generation
# ==================================
rule convert_sparse_matrix_to_dense:
    input:
        sparse_matrix="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/counts_sparse_matrix.txt",
        cells="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/cells.txt",
        genes="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/genes.txt"
    output:
        dense_matrix="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/counts_dense_matrix.txt",
    params:
        cores="1",
        hvmem="h_vmem=200G",
        memrequested="mem_requested=200G",
        tmprequested="tmp_requested=10G",
        tmpfree="tmpfree=10G",
        jobname="convert_sparse_matrix_to_dense"
    script:
        "code/celltype_deconvolution/convert_sparse_matrix_to_dense.snakemake.py"

# ==================================
# generate_pseudobulk_mixture_file
# ----------------------------------
# Generate pseudo_bulk mixture file from scRNA-Seq
# ==================================
rule generate_pseudobulk_mixture_file:
    input:
        seurat=lambda wildcards: expand("{analysis_data}", analysis_data=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020[wildcards.analysis])
    output:
        sum_matching_pseudobulk_mixture_matrix="analysis/CIBERSORTx/{analysis}/pseudobulk_mixture_matrix/sum_matching_pseudobulk_mixture_matrix.txt",
        avg_matching_pseudobulk_mixture_matrix="analysis/CIBERSORTx/{analysis}/pseudobulk_mixture_matrix/avg_matching_pseudobulk_mixture_matrix.txt"
    params:
        cores="1",
        hvmem="h_vmem=200G",
        memrequested="mem_requested=200G",
        tmprequested="tmp_requested=10G",
        tmpfree="tmpfree=10G",
        jobname="generate_pseudobulk_mixture_file"
    script:
        "code/celltype_deconvolution/generate_pseudobulk_mixture_file.snakemake.R"
