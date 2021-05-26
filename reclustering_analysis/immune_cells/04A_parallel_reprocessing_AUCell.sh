# Title: parallel reclustering script for each major lineage in the breast cancer mini-atlas dataset
#
# author: Sunny Wu
# group: TumourProgression
# version: 2.0.0
# date last edited: 20191013
#
# SUMMARY
#
# GLOBS
## Activate conda R enviroment
  source activate Renv
## R path
  Rdev="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R"
## Celltype IDs for reclustering
 CELLTYPES="T_cells Myeloid_cells B_cells Plasma_cells"
#  CELLTYPES="B_cells"
# STEPS TO RUN
    # REPROCESSING
    STEP01="TRUE"
# DIRECTORIES
## working directory
  TEMPPWD=$(pwd)

# 01 SUBMIT SUBSET REPROCESSING JOBS
if [[ "${STEP01}" == "TRUE" ]]; then
  echo "STEP01: Reprocessing Jobs"
  for celltype in ${CELLTYPES}; do
        JOBNAME="J04_${celltype}"
        CORES="8"
        MEMPERCORE="10G"
        if [[ "${celltype}" == "T_cells" ]]; then
          CORES="12"
          MEMPERCORE="10G"
          fi
        TEMPCELLTYPEDIR="${TEMPPWD}/analysis_01_${celltype}/"

      qsub \
      -cwd \
      -pe smp ${CORES} \
      -l mem_requested=${MEMPERCORE} \
      -b y \
      -j y \
      -N ${JOBNAME} \
      -V \
      -P TumourProgression \
      "${Rdev} CMD BATCH \
      --no-save '--args \
      ${celltype} \
      ${CORES} \
      ${MEMPERCORE} \
      ${TEMPCELLTYPEDIR}' \
      ${TEMPPWD}/08_reprocessing_script.R" \
      ${TEMPPWD}/analysis_01_${celltype}/04_reprocessing_log_${celltype}.txt \
      > ${TEMPPWD}/analysis_01_${celltype}/JOBID_04_reprocessing_${celltype}.txt

      echo "  submitted job ${JOBNAME}"
      done
    fi
#
