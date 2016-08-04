#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/SS_%j.out
#SBATCH --output=slurm/SS_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

SAMPLE=${1}
 INPUT=${2}
OUTPUT=${3}

HEADER="SS"

echo "$HEADER: ${SAMPLE} ${INPUT} to ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile; then exit 1; fi
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

module load ${MOD_JAVA}

PIC_ARGS="SORT_ORDER=coordinate \
CREATE_INDEX=true \
COMPRESSION_LEVEL=9 \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
TMP_DIR=${JOB_TEMP_DIR}"

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} SortSam ${PIC_ARGS} INPUT=${INPUT} OUTPUT=${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a ../commands.txt

${CMD}
if cmdFailed; then exit 1; fi

# Move output to final location
if ! finalOut; then exit 1; fi

rm ${INPUT} && echo "$HEADER: Purged aligned files!"

touch ${OUTPUT}.done

storeMetrics run
