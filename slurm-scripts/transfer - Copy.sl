#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		TransferFile
#SBATCH --time			0-00:20:00
#SBATCH --mem-per-cpu	512
#SBATCH --cpus-per-task	1
#SBATCH --error			slurm/TF_%j.out
#SBATCH --output		slurm/TF_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

   IDN=${1}
 INPUT=${2}
OUTPUT=${3}

HEADER="TF"

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
	echo "$HEADER: $LABEL already started!"
	exit 0
fi

if [ "${OUTPUT}" == "$(basename ${OUTPUT})" ]
then
	echo "$HEADER: Output file has not pathing. Storing in temp folder."
	OUTPUT="~/projects/NeSI_Transfer/${IDN}/${OUTPUT}"
fi

echo "$HEADER: ${INPUT} -> ${OUTPUT}"
date

CMD="ssh globus transfer --encrypt --perf-cc 4 --perf-p 8 --label \"${LABEL}\" -- nz#uoa/$(pwd)/${INPUT} nesi#otago-dtn01/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit $EXIT_PR
fi

touch ${LABEL}.transfer.done
