#!/bin/bash
#SBATCH --job-name		BaseRecalibrator
#SBATCH --time			0-02:00:00
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	8
#SBATCH --constraint	avx
#SBATCH --array			1-84
#SBATCH --error			slurm/BR_%A_%a.out
#SBATCH --output		slurm/BR_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

CONTIG=${CONTIGBLOCKS[$SLURM_ARRAY_TASK_ID]}
 INPUT=markdup/${CONTIG}.bam
OUTPUT=baserecal/${CONTIG}.firstpass

HEADER="BR"

echo "$HEADER: ${INPUT} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

GATK_PROC=BaseRecalibrator
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-knownSites ${DBSNP} \
-knownSites ${MILLS} \
-L ${CONTIG} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="srun $(which java) ${JAVA_ARGS} -jar ${GATK} ${GATK_ARGS} -I ${INPUT} -o ${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

JOBSTEP=0

if ! ${CMD}; then
	cmdFailed $?
	exit ${JOBSTEP}${EXIT_PR}
fi

storeMetrics

touch ${OUTPUT}.done
