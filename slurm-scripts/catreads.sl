#!/bin/bash
#SBATCH --job-name		CatReads
#SBATCH --time			0-01:00:00
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	1
#SBATCH --constraint	avx
#SBATCH --error			slurm/CR_%j.out
#SBATCH --output		slurm/CR_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

FILES=${1}
OUTPUT=${2}
IDN=$(echo $SLURM_JOB_NAME | cut -d'_' -f2)

BAMHEAD=$(echo ${FILES} | awk '{print $1}')

HEADER="CR"

echo $HEADER: $FILES + "Header($BAMHEAD) ->" $OUTPUT
date

for INPUT in ${FILES}; do
	if ! inFile; then exit $EXIT_IO; fi
done

if ! outDirs; then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

CMD="$(which srun) ${SAMTOOLS} cat -h ${BAMHEAD} -o ${JOB_TEMP_DIR}/${OUTPUT} ${FILES}"
echo "$HEADER: ${CMD}}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit $EXIT_PR
fi

# Move output to final location
if ! finalOut; then exit $EXIT_MV; fi

touch ${OUTPUT}.done

sbatch $(dispatch "RI") -J RI_${IDN} ${SLSBIN}/catreadsindex.sl ${OUTPUT} ${OUTPUT%.bam}.bai

if ! ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}; then
	echo "$HEADER: Transfer failed!"
	exit $EXIT_TF
fi
