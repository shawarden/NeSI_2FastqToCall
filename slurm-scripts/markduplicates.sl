#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=8192
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/MD_%A_%a.out
#SBATCH --output=slurm/MD_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
 INPUT=merged/${CONTIG}.bam
OUTPUT=markdup/${CONTIG}.bam

HEADER="MD"

echo "$HEADER: ${INPUT} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit 1; fi
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

PIC_ARGS="CREATE_INDEX=true \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
COMPRESSION_LEVEL=9 \
METRICS_FILE=${OUTPUT%.bam}.metrics \
TMP_DIR=${JOB_TEMP_DIR}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} MarkDuplicates ${PIC_ARGS} INPUT=${INPUT} OUTPUT=${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 1
fi

# Move output to final location
if ! finalOut; then exit 1; fi

rm ${INPUT} ${INPUT%.bam}.bai && echo "$HEADER: Purged input files!"

touch ${OUTPUT}.done

storeMetrics
