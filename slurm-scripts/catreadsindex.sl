#!/bin/bash
#SBATCH --job-name		ReadIndex
#SBATCH --time			0-01:30:00
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	1
#SBATCH --constraint	avx
#SBATCH --error			slurm/RI_%j.out
#SBATCH --output		slurm/RI_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

INPUT=${1}
OUTPUT=$([ "${2}" == "" ] && echo -ne "${INPUT%.bam}.bai" || echo -ne "${2}")
IDN=$(echo $SLURM_JOB_NAME | cut -d'_' -f2)

HEADER=RI

echo "$HEADER: ${INPUT} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit $EXIT_IO; fi
if ! outDirs; then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

CMD="$(which srun) ${SAMTOOLS} index ${INPUT} ${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit $EXIT_PR
fi

# Move output to final location
if ! finalOut; then exit $EXIT_MV; fi

touch ${OUTPUT}.done

if ! ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}; then
	echo "$HEADER: Transfer index failed!"
	exit $EXIT_TF
fi
