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
  # R34="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.4.1/bin/R"
## Celltype IDs for reclustering
  CELLTYPES="T_cells Myeloid_cells B_cells Plasma_cells"
# STEPS TO RUN
    # REPROCESSING
    STEP01="TRUE"
# DIRECTORIES
## working directory
  TEMPPWD=$(pwd)
## output directories
if [[ "${STEP01}" == "TRUE" ]]; then
  cd ${TEMPPWD}
  for celltype in ${CELLTYPES}; do
    mkdir "analysis_01_${celltype}"
    done
fi
#
# 01 SUBMIT SUBSET REPROCESSING JOBS
if [[ "${STEP01}" == "TRUE" ]]; then
  echo "STEP01: Reprocessing Jobs"
  for celltype in ${CELLTYPES}; do
        JOBNAME="J01_${celltype}"
        CORES="8"
        MEMPERCORE="20G"
        if [[ "${celltype}" == "T_cells" ]]; then
          CORES="24"
          MEMPERCORE="20G"
          fi
        TEMPCELLTYPEDIR="${TEMPPWD}/analysis_01_${celltype}"

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
      $TEMPCELLTYPEDIR' \
      ${TEMPPWD}/02_reprocessing_script.R" \
      ${TEMPPWD}/analysis_01_${celltype}/01_reprocessing_log_${celltype}.txt \
      > ${TEMPPWD}/analysis_01_${celltype}/JOBID_01_reprocessing_${celltype}.txt

      echo "  submitted job ${JOBNAME}"
      done
    fi
#
