# ==================================
# Dan Roden
# CREATED: 01/12/2020
# MODIFIED: 01/12/2020
# ----------------------------------
# Snakemake file for running DWLS decon
#
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

CELL_TYPES=["celltype_major", "celltype_minor", "celltype_subset"]

CANCER_CELL_ANNOTATIONS=["SCTyper_Dec"]

CYCLING_CELL_ANNOTATIONS=["combined"]

BULK_DATASETS={"pseudo_bulk/sum_matching": "analysis/CIBERSORTx/Jul2020/pseudobulk_mixture_matrix/sum_matching_pseudobulk_mixture_matrix.txt", "bulk_FFPE": "data/bulk_FFPE/counts/Aus-RiboZero_gene_rawcounts.matrix.txt", "metabric/discovery": "data/METABRIC/heloisa/METABRIC_Discovery.with_gene_heading.txt", "metabric/validation": "data/METABRIC/heloisa/METABRIC_Validation.with_gene_heading.txt"}

# ==================================
rule all:
    input:
        # dwls_create_signature_matrix
        expand("analysis/DWLS/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/signature_matrix/signature_matrix.txt", analysis=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020.keys(), normal_epithelial_cells=NORMAL_EPITHELIAL_CELL_CALLS, cell_type=CELL_TYPES, cancer_cell_annotation=CANCER_CELL_ANNOTATIONS, cycling_cell_annotation=CYCLING_CELL_ANNOTATIONS, subsample_fraction=SUBSAMPLE_FRACTION),
        # generate_pseudobulk_mixture_file
        expand("analysis/DWLS/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/cell_fractions/{bulk_dataset}/celltype_fractions.txt", analysis=ANALYSIS_ID_SEURAT_OBJ_MAPPINGS_2020.keys(), normal_epithelial_cells=NORMAL_EPITHELIAL_CELL_CALLS, cell_type=CELL_TYPES, cancer_cell_annotation=CANCER_CELL_ANNOTATIONS, cycling_cell_annotation=CYCLING_CELL_ANNOTATIONS, subsample_fraction=SUBSAMPLE_FRACTION, bulk_dataset=BULK_DATASETS.keys())

# ==================================
# dwls_create_signature_matrix
# ----------------------------------
# Create cell-type signature matrix (using MAST) for DWLS cell-type deconvolution
# ==================================
rule dwls_create_signature_matrix:
    input:
        count_matrix="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/counts_dense_matrix.txt",
        cell_ids="analysis/CIBERSORTx/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/cells.txt"
    output:
        signature_matrix="analysis/DWLS/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/signature_matrix/signature_matrix.txt"
    params:
        cores="1",
        hvmem="h_vmem=300G",
        memrequested="mem_requested=100G",
        tmprequested="tmp_requested=10G",
        tmpfree="tmpfree=10G",
        jobname="dwls_create_signature_matrix"
    log:
        "logs/ecotypes/DWLS/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/signature_matrix/output.log"
    script:
        "code/celltype_deconvolution/dwls_create_signature_matrix.snakemake.R"

# ==================================
# dwls_deconvolute_bulk
# ----------------------------------
# Run DWLS cell-type deconvolution on bulk samples
# ==================================
rule dwls_deconvolute_bulk:
    input:
        signature_matrix="analysis/DWLS/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/signature_matrix/signature_matrix.txt",
        bulk_data=lambda wildcards: expand("{bulk_dataset}", bulk_dataset=BULK_DATASETS[wildcards.bulk_dataset])
    output:
        celltype_fractions="analysis/DWLS/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/cell_fractions/{bulk_dataset}/celltype_fractions.txt"
    params:
        cores="1",
        hvmem="h_vmem=200G",
        memrequested="mem_requested=200G",
        tmprequested="tmp_requested=10G",
        tmpfree="tmpfree=10G",
        jobname="dwls_deconvolute_bulk"
    log:
        "logs/ecotypes/DWLS/{analysis}/{cell_type}/{normal_epithelial_cells}/{cancer_cell_annotation}/cycling_{cycling_cell_annotation}/sampled_{subsample_fraction}/cell_fractions/{bulk_dataset}/output.log"
    script:
        "code/celltype_deconvolution/dwls_deconvolute_bulk.snakemake.R"
