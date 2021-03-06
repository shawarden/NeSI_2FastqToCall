#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		Recalibration
#SBATCH --time			359
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	8
#SBATCH --constraint	avx
#SBATCH --array			1-84
#SBATCH --error			slurm/RC_%A_%a.out
#SBATCH --output		slurm/RC_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

CONTIG=$([ "$1" != "" ] && echo -ne "$1" || echo "${CONTIGBLOCKS[$SLURM_ARRAY_TASK_ID]}")
 INPUT=markdup/${CONTIG}.bam
  BQSR=bqsr_${CONTIG}.firstpass
OUTPUT=printreads/${CONTIG}.bam

HEADER="RC"

echo "$HEADER: ${INPUT} -> BQSR -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit $EXIT_IO; fi
if ! outDirs; then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

module load ${MOD_JAVA}
HEADER="BR"

CMD="srun $(which java) ${JAVA_ARGS} -jar ${GATK} ${GATK_BSQR} -L ${CONTIG} ${GATK_ARGS} -I ${INPUT} -o ${JOB_TEMP_DIR}/${BQSR}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

JOBSTEP=0
 
if ! ${CMD}; then
	cmdFailed $?
	exit ${JOBSTEP}${EXIT_PR}
fi

storeMetrics

BR_SECONDS=$SECONDS
SECONDS=0

HEADER="PR"

CMD="srun $(which java) ${JAVA_ARGS} -jar ${GATK} ${GATK_READ} -L ${CONTIG} ${GATK_ARGS} -I ${INPUT} -BQSR ${JOB_TEMP_DIR}/${BQSR} -o ${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

JOBSTEP=1

if ! ${CMD}; then
	cmdFailed $?
	exit ${JOBSTEP}${EXIT_PR}
fi

storeMetrics

SECONDS=$(($SECONDS + $BR_SECONDS))
JOBSTEP=""

rm ${INPUT} ${INPUT%.bam}.bai && echo "$HEADER: Purged input files!"

touch ${OUTPUT}.done
