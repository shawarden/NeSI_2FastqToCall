#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		MergeContigs
#SBATCH --time			0-02:00:00
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	4
#SBATCH --constraint	avx
#SBATCH --array			1-84
#SBATCH --error			slurm/MC_%A_%a.out
#SBATCH --output		slurm/MC_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

CONTIG=${CONTIGARRAY[$SLURM_ARRAY_TASK_ID]}
OUTPUT=merged/${CONTIG}.bam

HEADER="MC"

echo "$HEADER: ${CONTIG} -> ${OUTPUT}"
date

contigMerBlocks=$(find . -type f -iwholename "*/split/${CONTIG}_*.bam" -printf '%h\0%d\0%p\n' | sort -t '\0' -n | awk -F '\0' '{print $3}')
numcontigMerBlocks=$(echo "$contigMerBlocks" | wc -l)

if [ $numcontigMerBlocks -eq 0 ]; then
	echo "$HEADER: Merge contig ${CONTIG} contains $numcontigMerBlocks files!"
#	scriptFailed
	exit $EXIT_IO
else
	echo $HEADER: Merge contig ${CONTIG} will run $numcontigMerBlocks files: \"${contigMerBlocks}\"
fi

mergeList=""
for INPUT in ${contigMerBlocks}; do
	if inFile; then
		mergeList="${mergeList} I=${INPUT}"
	else
		exit $EXIT_IO
	fi
done

if [ "$mergeList" == "" ]; then
	echo "$HEADER: No inputs defined!"
	exit $EXIT_IO
fi

# Make sure input and target folders exists and that output file does not!
if ! outDirs; then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

PIC_ARGS="SORT_ORDER=coordinate \
USE_THREADING=true \
CREATE_INDEX=true \
COMPRESSION_LEVEL=9 \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
TMP_DIR=${JOB_TEMP_DIR}"

module load ${MOD_JAVA}

CMD="srun $(which java) ${JAVA_ARGS} -jar ${PICARD} MergeSamFiles ${PIC_ARGS} ${mergeList} O=${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

JOBSTEP=0

if ! ${CMD}; then
	cmdFailed $?
	exit ${JOBSTEP}${EXIT_PR}
fi

storeMetrics

# Move output to final location
if ! finalOut; then exit $EXIT_MV; fi

# Get list of sorted blocks and add them to the delete pile!
# Don't delete the .done files through.
sortedBlocks=$(find . -type f -iwholename "*/sorted/*.ba*" | grep -v ".done")
numSortedBlocks=$(echo "$sortedBlocks" | wc -l)

if [ $numSortedBlocks -gt 0 ]; then
	rm ${sortedBlocks} && echo "$HEADER: Purged $numSortedBlocks sorted block files!"
fi

rm ${contigMerBlocks} && echo "$HEADER: Purging $numcontigMerBlocks contig merge blocks!"

touch ${OUTPUT}.done
