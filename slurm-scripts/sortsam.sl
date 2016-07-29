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

echo "SS: ${SAMPLE} ${INPUT} to ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "SS: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "SS"
	exit 1
fi

module load ${MOD_JAVA}

PIC_ARGS="SORT_ORDER=coordinate \
CREATE_INDEX=true \
COMPRESSION_LEVEL=9 \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
TMP_DIR=${JOB_TEMP_DIR}"

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} SortSam ${PIC_ARGS} INPUT=${INPUT} OUTPUT=${OUTPUT}"
echo "SS: ${CMD}" | tee -a ../commands.txt

${CMD}
passed=$?

echo "SS: Sorted ${INPUT} to ${OUTPUT} in $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "SS: Failed ${SAMPLE} ${INPUT}!"
#	scriptFailed "SS"
	exit 1
fi

rm ${INPUT}

touch ${OUTPUT}.done

storeMetrics run
