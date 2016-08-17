#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		SortSAM
#SBATCH --time			0-01:00:00
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	4
#SBATCH --constraint	avx
#SBATCH --array			0-999
#SBATCH --error			slurm/SS_%A_%a.out
#SBATCH --output		slurm/SS_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

SAMPLE=${1}

 BLOCK=$(printf "%0${FASTQ_MAXZPAD}d" $SLURM_ARRAY_TASK_ID)
 INPUT=aligned/${BLOCK}.bam
OUTPUT=sorted/${BLOCK}.bam

HEADER="SS"

echo "$HEADER: ${SAMPLE} ${INPUT} to ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit $EXIT_IO; fi
if ! outDirs; then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

module load ${MOD_JAVA}

PIC_ARGS="SORT_ORDER=coordinate \
CREATE_INDEX=true \
COMPRESSION_LEVEL=9 \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
TMP_DIR=${JOB_TEMP_DIR}"

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} SortSam ${PIC_ARGS} INPUT=${INPUT} OUTPUT=${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a ../commands.txt

if ! ${CMD}; then
	cmdFailed
	exit $EXIT_PR
fi

# Move output to final location
if ! finalOut; then exit $EXIT_MV; fi

rm ${INPUT} && echo "$HEADER: Purged aligned files!"

touch ${OUTPUT}.done

storeMetrics run
