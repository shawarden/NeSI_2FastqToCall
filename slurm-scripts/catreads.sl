#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/CR_%j.out
#SBATCH --output=slurm/CR_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

FILES=${1}
OUTPUT=${2}
IDN=$(echo $SLURM_JOB_NAME | cut -d'_' -f2)

BAMHEAD=$(echo ${INPUT} | awk '{print $1}')

HEADER="CR"

echo $HEADER: $FILES + "Header($BAMHEAD) ->" $OUTPUT
date

for INPUT in ${FILES}; do
	if ! inFile; then exit 1; fi
done

if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

CMD="$(which srun) ${SAMTOOLS} cat -h ${BAMHEAD} -o ${JOB_TEMP_DIR}/${OUTPUT} ${FILES}"
echo "$HEADER: ${CMD1}, ${CMD2}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 1
fi

# Move output to final location
if ! finalOut; then exit 1; fi

touch ${OUTPUT}.done

sbatch -J TR_${IDN} ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}
sbatch -J RI_${IDN} ${SLSBIN}/catreadsindex.sl ${OUTPUT} ${OUTPUT%.bam}.bai

storeMetrics


