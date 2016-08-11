#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		PrintReads
#SBATCH --time			0-03:00:00
#SBATCH --mem-per-cpu	4092
#SBATCH --cpus-per-task	8
#SBATCH --constraint	avx
#SBATCH --array			1-84
#SBATCH --error			slurm/PR_%A_%a.out
#SBATCH --output		slurm/PR_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
 INPUT=markdup/${CONTIG}.bam
  BQSR=baserecal/${CONTIG}.firstpass
OUTPUT=printreads/${CONTIG}.bam

HEADER="PR"

echo "$HEADER: ${INPUT} + ${BQSR} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit 1; 
if ! (INPUT=$(echo $BQSR); inFile); then exit 1
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

GATK_PROC=PrintReads
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
--bam_compression 9 \
-L ${CONTIG} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${GATK} ${GATK_ARGS} -I ${INPUT} -BQSR ${BQSR} -o ${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 1
fi

# Move output to final location
if ! finalOut; then exit 1; fi

rm ${INPUT} ${INPUT%.bam}.bai ${INPUT%.bam}.metrics && echo "$HEADER: Purged input files!"

touch ${OUTPUT}.done

storeMetrics
