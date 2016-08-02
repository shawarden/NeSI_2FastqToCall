#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=8192
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/MD_%j.out
#SBATCH --output=slurm/MD_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
 INPUT=merged/${CONTIG}.bam
OUTPUT=markdup/${CONTIG}.bam

echo "MD: ${INPUT} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "MD: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "MD"
	exit 1
fi

PIC_ARGS="CREATE_INDEX=true \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
COMPRESSION_LEVEL=9 \
METRICS_FILE=${OUTPUT%.bam}.metrics \
TMP_DIR=${JOB_TEMP_DIR}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} MarkDuplicates ${PIC_ARGS} INPUT=${INPUT} OUTPUT=${OUTPUT}"
echo "MD: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "MD: ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]
then
	echo "MD: ${INPUT} failed!"
#	scriptFailed "MD"
	exit 1
fi

# Remove old data
rm ${INPUT} ${INPUT%.bam}.bai

touch ${OUTPUT}.done

storeMetrics
