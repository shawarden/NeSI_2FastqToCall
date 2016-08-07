#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --array=0-999
#SBATCH --error=slurm/PA_%A_%a.out
#SBATCH --output=slurm/PA_%A_%a.out

source /projects/uoo00032/Resources/bin/baserefs.sh

   SAMPLE=${1}
READGROUP=${2}
    BLOCK=$(printf "%0${FASTQ_MAXZPAD}d" $SLURM_ARRAY_TASK_ID)
   
   OUTPUT=aligned/${BLOCK}.bam

HEADER="PA"
   
READ1=blocks/${READGROUP}_R1_${BLOCK}.fastq.gz
READ2=blocks/${READGROUP}_R2_${BLOCK}.fastq.gz

echo "$HEADER: $READGROUP $READ1 $READ2 -> $OUTPUT"
date

# Make sure input and target folders exists and that output file does not!
if ! (INPUT=${READ1}; inFile); then exit 1
if ! (INPUT=${READ2}; inFile); then exit 1
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

# Get readgroup blocks from either INFO_INFO_INFO_.. or INFO INFO INFO ...
     INTRUMENT=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $1}')
INSTRUMENT_RUN=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $2}')
     FLOW_CELL=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $3}')
     CELL_LANE=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $4}')
         INDEX=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $5}')

RG_ID="ID:${INTRUMENT}_${INSTRUMENT_RUN}_${FLOW_CELL}_${CELL_LANE}_${INDEX}"
RG_PL="PL:Illumina"
RG_PU="PU:${FLOW_CELL}.${CELL_LANE}"
RG_LB="LB:${SAMPLE}"
RG_SM="SM:$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')"

echo "$HEADER: ${BWA} mem -M -t ${SLURM_JOB_CPUS_PER_NODE} -R @RG'\t'$RG_ID'\t'$RG_PL'\t'$RG_PU'\t'$RG_LB'\t'$RG_SM $REF $READ1 $READ2 | ${SAMTOOLS} view -bh - > ${OUTPUT}" | tee -a ../commands.txt

${BWA} mem -M -t ${SLURM_JOB_CPUS_PER_NODE} -R @RG'\t'$RG_ID'\t'$RG_PL'\t'$RG_PU'\t'$RG_LB'\t'$RG_SM $REF $READ1 $READ2 | ${SAMTOOLS} view -bh - > ${JOB_TEMP_DIR}/${OUTPUT}
if cmdFailed; then exit 1; fi

# Move output to final location
if ! finalOut; then exit 1; fi

rm ${READ1} ${READ2} && echo "$HEADER: Purged read files!"

touch ${OUTPUT}.done

storeMetrics run
