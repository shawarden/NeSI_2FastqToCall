#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		TransferFile
#SBATCH --time			359
#SBATCH --mem-per-cpu	512
#SBATCH --cpus-per-task	1
#SBATCH --error			slurm/TF_%j.out
#SBATCH --output		slurm/TF_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

   IDN=${1}
 INPUT=${2}
OUTPUT=${3}

echo "$HEADER: IDN:${IDN} + INPUT:${INPUT} -> OUTPUT:${OUTPUT}"

if ! inFile; then exit $EXIT_IO; fi

if [ "${OUTPUT}" == "" ]
then
	echo "$HEADER: Output not defined. Using input:${INPUT}"
	OUTPUT=${INPUT}
fi

LABEL=$(echo $(basename ${OUTPUT}) | tr '.' '_')

if [ -e ${LABEL}.transfer.done ]
then
	echo "$HEADER: $LABEL already done!"
	exit 0
fi

# Strip pathing from output
OUTPUT="~/projects/SequenceData/GRCh37/Alignments/Genomic/${IDN}/$(basename ${OUTPUT})"

echo "$HEADER: ${INPUT} -> ${OUTPUT}"
date

CMD="globus transfer --label \"${LABEL}\" -- ${ENDPOINT_NESI}:/$(pwd)/${INPUT} ${ENDPOINT_UOO}:/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

JOBSTEP=""

if ! globus task wait $(${CMD} | grep "Task ID" | awk '{print $3}'); then
	cmdFailed $?
	exit $EXIT_PR
fi

storeMetrics

touch ${LABEL}.transfer.done
