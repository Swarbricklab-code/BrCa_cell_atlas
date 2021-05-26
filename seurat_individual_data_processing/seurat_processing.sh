# SEURAT V3 JOB SUBMISSION SCRIPT
# this script is for the processing of individual seurat objects
# SUNNY WU
#
#
# Activate conda R enviroment
source activate Renv
#
# R PATH
R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R"
TEMPPWD=$(pwd)
mkdir output/
# make output directories for individual samples
for samplename in $(cut -d ',' -f 2 ./config/input_sample_file.csv | tail -n +2); do
  mkdir output/seurat_$samplename/
done
# submit jobs to cluster
for samplenum in $(cut -d ',' -f 1 $TEMPPWD/config/input_sample_file.csv | tail -n +2); do

  SAMPLENAME=$(cut -d ',' -f 2 $TEMPPWD/config/input_sample_file.csv | tail -n +2 | head -n $samplenum | tail -n 1)
  REFERENCEGENOME="human"
  SEURATGENEINPUTFILE="${TEMPPWD}/config/seurat_gene_input_file.csv"
  SEURATPARAMFILE="${TEMPPWD}/config/seurat_params_file.csv"
  SEURATSCRIPT="${TEMPPWD}/config/scripts/seurat_HPC_processing.R"
  SEURATJOBNAME="s_${SAMPLENAME}"
  MATRIXPATH="/paella/TumourProgressionGroupTemp/projects/data/cellranger_count/human_breast/${SAMPLENAME}/count_*GRCh38/outs/raw_gene_bc_matrices/GRCh38/"
  TEMPWD="${TEMPPWD}/output/seurat_${SAMPLENAME}/"

  cd ${TEMPPWD}/output/seurat_${SAMPLENAME}/

  qsub \
  -cwd \
  -pe smp 8 \
  -l mem_requested=10G \
  -b y \
  -j y \
  -V \
  -N ${SEURATJOBNAME}\
  -P TumourProgression \
  "${R} CMD BATCH \
  --no-save '--args \
  ${SAMPLENAME} \
  ${REFERENCEGENOME} \
  ${TEMPWD} \
  ${MATRIXPATH} \
  ${SEURATGENEINPUTFILE} \
  ${SEURATPARAMFILE}' \
  ${SEURATSCRIPT}"

  cd ${TEMPPWD}

done
