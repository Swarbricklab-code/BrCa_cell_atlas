#!/bin/bash
#
# 10X Cellranger pipeline
# script name: sunwu_10X_pipeline.sh
# written by; Sunny Wu
# group; Cancer Tumour Progression Group Garvan Institute
#
# specify input files in config/ folder
#     - config/sample_input_file.csv
#               path_to_bclfiles,/share/ScratchGeneral/sunwu/bclfiles/Files/,,,,
#               sample_num,sampleID,barcodes,ref_species,cell_number
#               1,CID4404,SI-GA-A1,human,9000
#               2,CID4398,SI-GA-A2,human,9000
#      - config/cellranger_input_samplesheet.csv
#                cellranger mkfastq input template file
#
# IMPORTANT: Execute this script in a screen using;
#   screen -S projectname
#   sh sunwu_10X_pipeline.sh
#
# STEPS TO RUN
    # MKDIR
    STEP01="TRUE"
    # MKFASTQ
    STEP02="TRUE"
    # CELLRANGER COUNT
    STEP03="TRUE"
#
##--------------------------------------- GLOBAL VARIABLES
### Path to BCL files
BCLPATH=$(cut -d ',' -f 2 ./config/sample_input_file.csv | head -n 1)
echo "bcl path"
echo ${BCLPATH}
###  working directory
TEMPPWD=$(pwd)
###  path to cellranger v2.2.0
CELLRANGER="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/cellranger-3.0.2/cellranger"
##--------------------------------------- MAKE DIRECTORIES
if [[ "${STEP01}" == "TRUE" ]]; then

  mkdir output
  mkdir output/mkfastq/
  mkdir output/count/
  mkdir output/seurat/

  ### sub directories for each sample
  for samplename in $(cut -d ',' -f 2 ./config/sample_input_file.csv | tail -n +3); do
    echo "making directory for $samplename"
    mkdir output/count/$samplename/
    mkdir output/seurat/seurat_$samplename/
  done
fi
##--------------------------------------- CELLRANGER MKFASTQ
if [[ "${STEP02}" == "TRUE" ]]; then
    ### generate input spreadsheet
    for samplenum in $(cut -d ',' -f 1 ./config/sample_input_file.csv | tail -n +3); do
      SAMPLENAME=$(cut -d ',' -f 2 ./config/sample_input_file.csv | tail -n +3 | head -n $samplenum | tail -n 1)
      SAMPLEBARCODE=$(cut -d ',' -f 3 ./config/sample_input_file.csv | tail -n +3 | head -n $samplenum | tail -n 1)
      CSVSTRING="1,$SAMPLENAME,$SAMPLENAME,,,$SAMPLEBARCODE,$SAMPLEBARCODE,Clinical_10X,"
      CSVSTRING2="2,$SAMPLENAME,$SAMPLENAME,,,$SAMPLEBARCODE,$SAMPLEBARCODE,Clinical_10X,"
      CSVSTRING3="3,$SAMPLENAME,$SAMPLENAME,,,$SAMPLEBARCODE,$SAMPLEBARCODE,Clinical_10X,"
      CSVSTRING4="4,$SAMPLENAME,$SAMPLENAME,,,$SAMPLEBARCODE,$SAMPLEBARCODE,Clinical_10X,"
      CSVSTRINGFINAL="${CSVSTRING}\n${CSVSTRING2}\n${CSVSTRING3}\n${CSVSTRING4}"
      echo -e "${CSVSTRINGFINAL}"
    done > temp_file.csv
    cat ./config/cellranger_input_samplesheet.csv ./temp_file.csv > output/mkfastq/cellranger_input_samplesheet.csv
    rm temp_file.csv

    cd output/mkfastq/
    ### submit mkfastq job
      qsub \
      -cwd \
      -pe smp 16 \
      -l mem_requested=20G \
      -P TumourProgression \
      -b y \
      -j y \
      -V \
      -N mkfastq \
      $CELLRANGER mkfastq \
      --run=$BCLPATH \
      --samplesheet=./cellranger_input_samplesheet.csv \
      --lanes=1,2,3,4 \
      --id=mkfastq_output \
      --qc \
      > ./mkfastq_job_id.txt

    ### job ID for subsequent qsub holds
    MKFASTQJOBID=$(cut -d ' ' -f 3 ./mkfastq_job_id.txt)

    echo "     ++++ cellranger mkfastq job submission"
    cat ./mkfastq_job_id.txt
fi
##--------------------------------------- CELLRANGER COUNT
if [[ "${STEP03}" == "TRUE" ]]; then
    echo "     ++++ cellranger count pipeline"
    for samplenum in $(cut -d ',' -f 1 ${TEMPPWD}/config/sample_input_file.csv | tail -n +3); do
      cd $TEMPPWD

      SAMPLENAME=$(cut -d ',' -f 2 ${TEMPPWD}/config/sample_input_file.csv | tail -n +3 | head -n $samplenum | tail -n 1)
      SAMPLEBARCODE=$(cut -d ',' -f 3 ${TEMPPWD}/config/sample_input_file.csv | tail -n +3 | head -n $samplenum | tail -n 1)
      CELLNUMBER=$(cut -d ',' -f 5 ${TEMPPWD}/config/sample_input_file.csv | tail -n +3 | head -n $samplenum | tail -n 1)
      REFERENCEGENOME=$(cut -d ',' -f 4 ${TEMPPWD}/config/sample_input_file.csv | tail -n +3 | head -n $samplenum | tail -n 1)
      FASTQPATH="$TEMPPWD/output/mkfastq/mkfastq_output/outs/fastq_path/Clinical_10X/$SAMPLENAME/"
      COUNTJOBID="./countjobid_$SAMPLENAME.txt"
      COUNTJOBNAME="jidcount_$SAMPLENAME"

      # transcriptome reference path
      if [ "$REFERENCEGENOME" = "human" ]; then
        TRANSCRIPTOME="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/refdata/refdata-cellranger-GRCh38-1.2.0/"
      elif [ "$REFERENCEGENOME" = "mouse" ]; then
        TRANSCRIPTOME="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/refdata/refdata-cellranger-mm10-1.2.0/"
      elif [ "$REFERENCEGENOME" = "both" ]; then
        TRANSCRIPTOME="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/refdata/refdata-cellranger-GRCh38_mm10-1.2.0/"
      fi

      # genome reference output file name
      if [ "$REFERENCEGENOME" = "human" ]; then
        OUTPUTSTRINGREF="GRCh38"
      elif [ "$REFERENCEGENOME" = "mouse" ]; then
        OUTPUTSTRINGREF="mm10"
      elif [ "$REFERENCEGENOME" = "both" ]; then
        OUTPUTSTRINGREF="GRCh38_mm10"
      fi

      OUTPUT_ID_STRING="count_${SAMPLENAME}_${OUTPUTSTRINGREF}"

      # run cellranger count
      TEMPWD="$TEMPPWD/output/count/$SAMPLENAME/"
      cd $TEMPWD

      qsub \
      -cwd \
      -pe smp 16 \
      -l mem_requested=20G \
      -P TumourProgression \
      -b y \
      -j y \
      -V \
      -hold_jid $MKFASTQJOBID \
      -N $COUNTJOBNAME \
      "${CELLRANGER} count \
      --id=$OUTPUT_ID_STRING \
      --sample=$SAMPLENAME \
      --transcriptome=$TRANSCRIPTOME \
      --fastqs=$FASTQPATH \
      --indices=$SAMPLEBARCODE \
      --expect-cells=$CELLNUMBER" \
      > $COUNTJOBID

      echo "Sample:$SAMPLENAME"
      cat $COUNTJOBID

      cd $TEMPPWD
    done
fi
