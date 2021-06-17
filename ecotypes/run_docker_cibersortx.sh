#!/bin/bash
# ==========================
# Dan roden
# runs CIBERSORTx using docker
# ==========================
# ---------------
# Assuming input directory stucture of:
# ANALYSIS_DIR="analysis/CIBERSORTx/${ANALYSIS_ID}/${CELL_TYPE}/${NORMAL_EPITHELIAL_CELLS}/${CANCER_CELL_ANNOTATION}/cycling_${CYCLING_CELL_ANNOTATION}/sampled_${FRACTION_SAMPLED}/"
# ARGS:
# $1 = ANALYSIS_ID
# $2 = CELL_TYPE
# $3 = NORMAL_EPITHELIAL_CELLS
# $4 = CANCER_CELL_ANNOTATION
# $5 = CYCLING_CELL_ANNOTATION
# $6 = FRACTION_SAMPLED
# $7 = PERMUTATIONS
#
# ==========================
ROOT_DIR=""
echo ${ROOT_DIR}

# ==========================
# set params
# ==========================
ANALYSIS_ID=$1
echo "ANALYSIS_ID: "${ANALYSIS_ID}

CELL_TYPE=$2
echo "CELL_TYPE: "${CELL_TYPE}

NORMAL_EPITHELIAL_CELLS=$3
echo "NORMAL_EPITHELIAL_CELLS: "${NORMAL_EPITHELIAL_CELLS}

CANCER_CELL_ANNOTATION=$4
echo "CANCER_CELL_ANNOTATION: "${CANCER_CELL_ANNOTATION}

CYCLING_CELL_ANNOTATION=$5
echo "CYCLING_CELL_ANNOTATION: "${CYCLING_CELL_ANNOTATION}

FRACTION_SAMPLED=$6
echo "FRACTION_SAMPLED: "${FRACTION_SAMPLED}

PERMUTATIONS=$7
echo "PERMUTATIONS: "${PERMUTATIONS}

ANALYSIS_DIR="analysis/CIBERSORTx/${ANALYSIS_ID}/${CELL_TYPE}/${NORMAL_EPITHELIAL_CELLS}/${CANCER_CELL_ANNOTATION}/cycling_${CYCLING_CELL_ANNOTATION}/sampled_${FRACTION_SAMPLED}/"
echo "ANALYSIS_DIR: "${ANALYSIS_DIR}

# --------------------------
# Generate signature matrix
# --------------------------
mkdir -p "${ANALYSIS_DIR}/signature_matrix/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/"

sudo docker run --name ${CELL_TYPE}_sigmatrix -v "${ROOT_DIR}${ANALYSIS_DIR}":/src/data -v "${ROOT_DIR}${ANALYSIS_DIR}signature_matrix/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/":/src/outdir cibersortx/fractions --username d.roden@garvan.org.au --token 55e2e1679b8027f30578d4320e1baf7a --refsample counts_dense_matrix.txt --single_cell TRUE --G.min 300 --G.max 500 --q.value 0.01 --filter FALSE --k.max 999 --remake TRUE --replicates 5 --sampling 0.5 --fraction 0.75 --verbose TRUE 2>&1 | tee "${ANALYSIS_DIR}signature_matrix/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/generate_signature_matrix.log"
sudo docker container rm ${CELL_TYPE}_sigmatrix

# --------------------------
# Calculate cell fractions
# --------------------------
# Copy needed files to temp "input staging" directory
# --------------------------------------------------------
# signature matrix
SIG_MATRIX="sigmatrix.${ANALYSIS_ID}.${CELL_TYPE}.${NORMAL_EPITHELIAL_CELLS}.${CANCER_CELL_ANNOTATION}.${CYCLING_CELL_ANNOTATION}.${FRACTION_SAMPLED}.txt"
cp "${ANALYSIS_DIR}signature_matrix/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/CIBERSORTx_counts_dense_matrix_inferred_phenoclasses.CIBERSORTx_counts_dense_matrix_inferred_refsample.bm.K999.txt" "analysis/CIBERSORTx/temp_input/${SIG_MATRIX}"

# reference matrix (sc counts)
REF_MATRIX="counts_dense_matrix.${ANALYSIS_ID}.${CELL_TYPE}.${NORMAL_EPITHELIAL_CELLS}.${CANCER_CELL_ANNOTATION}.${CYCLING_CELL_ANNOTATION}.${FRACTION_SAMPLED}.txt"
cp "${ANALYSIS_DIR}counts_dense_matrix.txt" "analysis/CIBERSORTx/temp_input/${REF_MATRIX}"

# --------------------------
# cell-fractions: pseudo-bulk
# --------------------------

# with no batch correction
mkdir -p "${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/pseudo_bulk/sum_matching/rmbatch_none/perm_${PERMUTATIONS}/relative/"

sudo docker run --name ${CELL_TYPE}_pseudobulk -v "${ROOT_DIR}analysis/CIBERSORTx/temp_input/":/src/data -v "${ROOT_DIR}${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/pseudo_bulk/sum_matching/rmbatch_none/perm_${PERMUTATIONS}/relative/":/src/outdir cibersortx/fractions --username d.roden@garvan.org.au --token 55e2e1679b8027f30578d4320e1baf7a --mixture "sum_matching_pseudobulk_mixture_matrix.txt" --sigmatrix ${SIG_MATRIX} --perm ${PERMUTATIONS} --rmbatchSmode FALSE --absolute FALSE 2>&1 | tee "${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/pseudo_bulk/sum_matching/rmbatch_none/perm_${PERMUTATIONS}/relative/cell_fractions.log"
sudo docker container rm ${CELL_TYPE}_pseudobulk


# --------------------------
# METABRIC
# --------------------------
# QN FALSE
QN="FALSE"

# --------------------------
# cell-fractions: METABRIC (discovery)
# --------------------------

# with S-mode batch correction
mkdir -p "${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/discovery/rmbatchSmode/refsample_counts/QN_${QN}/perm_${PERMUTATIONS}/relative/"

sudo docker run --name ${CELL_TYPE}_metadisc_smode -v "${ROOT_DIR}analysis/CIBERSORTx/temp_input/":/src/data -v "${ROOT_DIR}${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/discovery/rmbatchSmode/refsample_counts/QN_${QN}/perm_${PERMUTATIONS}/relative/":/src/outdir cibersortx/fractions --username d.roden@garvan.org.au --token 55e2e1679b8027f30578d4320e1baf7a --mixture "METABRIC_Discovery.with_gene_heading.txt" --sigmatrix ${SIG_MATRIX} --refsample ${REF_MATRIX} --perm ${PERMUTATIONS} --rmbatchSmode TRUE --absolute FALSE --QN ${QN} 2>&1 | tee "${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/discovery/rmbatchSmode/refsample_counts/QN_${QN}/perm_${PERMUTATIONS}/relative/cell_fractions.log"
sudo docker container rm ${CELL_TYPE}_metadisc_smode

# --------------------------
# cell-fractions: METABRIC (validation)
# --------------------------

# with S-mode batch correction
mkdir -p "${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/validation/rmbatchSmode/refsample_counts/QN_${QN}/perm_${PERMUTATIONS}/relative/"

sudo docker run --name ${CELL_TYPE}_metavalid_smode -v "${ROOT_DIR}analysis/CIBERSORTx/temp_input/":/src/data -v "${ROOT_DIR}${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/validation/rmbatchSmode/refsample_counts/QN_${QN}/perm_${PERMUTATIONS}/relative/":/src/outdir cibersortx/fractions --username d.roden@garvan.org.au --token 55e2e1679b8027f30578d4320e1baf7a --mixture "METABRIC_Validation.with_gene_heading.txt" --sigmatrix ${SIG_MATRIX} --refsample ${REF_MATRIX} --perm ${PERMUTATIONS} --rmbatchSmode TRUE --absolute FALSE --QN ${QN} 2>&1 | tee "${ANALYSIS_DIR}cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/validation/rmbatchSmode/refsample_counts/QN_${QN}/perm_${PERMUTATIONS}/relative/cell_fractions.log"
sudo docker container rm ${CELL_TYPE}_metavalid_smode
