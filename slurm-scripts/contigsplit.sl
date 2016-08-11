#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		ContigSplit
#SBATCH --time			0-00:10:00
#SBATCH --mem-per-cpu	2048
#SBATCH --cpus-per-task	2
#SBATCH --constraint	avx
#SBATCH --array			1-84
#SBATCH --error			slurm/CS_%A_%a.out
#SBATCH --output		slurm/CS_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

INPUT=${1}
BLOCK=${2}

HEADER="CS"

CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
OUTPUT=split/${CONTIG}_${BLOCK}.bam

echo "$HEADER: ${INPUT} by contig ${CONTIG} to ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit 1; fi
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

CMD="${SAMTOOLS} view -bh -@ 1 ${INPUT} ${CONTIG} > ${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a ../commands.txt

if ! eval ${CMD}; then
	cmdFailed
	exit 1
fi

# Move output to final location
if ! finalOut; then exit 1; fi

touch ${OUTPUT}.done

storeMetrics run
