#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		Packup
#SBATCH --time			0-01:00:00
#SBATCH --mem-per-cpu	2048
#SBATCH --cpus-per-task	1
#SBATCH --error			slurm/PU_%j.out
#SBATCH --output		slurm/PU_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

   IDN=${1}

HEADER="PU"

echo "$HEADER: IDN:${IDN}"

rm printreads/*.bam printreads/*.bai haplo/*.g.vcf.gz.tbi 

if ! ${SLSBIN}/transfer.sl ${IDN} metrics.txt ${IDN}.metrics.txt; then
	echo "$HEADER: Transfer failed!"
	exit $EXIT_TF
fi

if ! ${SLSBIN}/transfer.sl ${IDN} commands.txt ${IDN}.commands.txt; then
	echo "$HEADER: Transfer failed!"
	exit $EXIT_TF
fi



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
	cmdFailed $?
	exit $EXIT_PR
fi

touch ${LABEL}.transfer.done
