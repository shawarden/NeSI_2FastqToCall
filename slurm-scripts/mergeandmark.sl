#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		MergeAndMark
#SBATCH --time			359
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	8
#SBATCH --constraint	avx
#SBATCH --array			1-84
#SBATCH --error			slurm/MM_%A_%a.out
#SBATCH --output		slurm/MM_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

CONTIG=${CONTIGARRAY[$SLURM_ARRAY_TASK_ID]}
MERGED=${JOB_TEMP_DIR}/merged.bam
OUTPUT=markdup/${CONTIG}.bam

HEADER="MM"

echo "$HEADER: ${CONTIG} -> merged -> ${OUTPUT}"
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
		mergeList="${mergeList} INPUT=${INPUT}"
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

module load ${MOD_JAVA}

#CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} MergeSamFiles ${PIC_ARGS} ${MERGE_ARGS} ${mergeList} OUTPUT=${MERGED}"
#echo "$HEADER: ${CMD}" | tee -a commands.txt

#if ! ${CMD}; then
#	cmdFailed $?
#	exit ${EXIT_PR}
#fi

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} MarkDuplicates ${PIC_ARGS} ${MARK_ARGS} ${mergeList} OUTPUT=${OUTPUT}" #INPUT=${MERGED} 
echo "$HEADER: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed $?
	exit ${EXIT_PR}
fi

rm ${contigMerBlocks} && echo "$HEADER: Purging $numcontigMerBlocks contig merge blocks!"

touch ${OUTPUT}.done

