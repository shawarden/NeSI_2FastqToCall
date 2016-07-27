#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=4092
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/mergereads_%j.out
#SBATCH --output=slurm/mergereads_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

 INPUT=${1}
OUTPUT=${2}
 CLEAN=${3}

echo "MergeReads: ${INPUT} -> ${OUTPUT}"
date

inputCount=$(echo "${INPUT}" | wc -w)
contigCount=$(echo "${CONTIGS}" | wc -w)

if [ $inputCount -ne $contigCount ]; then
	echo "MergeReads: Input chunks (${inputCount}) doesn't match contig count (${contigCount})!"
#	scriptFailed "MergeReads"
	exit 1
else
	echo "MergeReads: Input chunks (${inputCount}) matches contig count (${contigCount})."
fi

mergeList=""
for file in ${INPUT}; do
	if [ -e $file ]; then
		mergeList="${mergeList} I=${file}"
	else
		echo "MergeReads: Input file $file doesn't exist!"
		scriptFailed "MergeReads"
		exit 1
	fi
done

if [ "$mergeList" == "" ]; then
	echo "MergeReads: No inputs defined!"
#	scriptFailed "MergeReads"
	exit 1
fi

PIC_ARGS="SORT_ORDER=coordinate \
USE_THREADING=true \
CREATE_INDEX=true \
COMPRESSION_LEVEL=9 \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
TMP_DIR=${JOB_TEMP_DIR}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} MergeSamFiles ${PIC_ARGS} ${mergeList} O=${OUTPUT}"
echo "MergeReads: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

# Report run time.
echo "MergeReads: ran for $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "MergeReads: ${CONTIG} failed!"
#	scriptFailed "MergeReads"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics
