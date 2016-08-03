#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=512
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --error=slurm/TF_%j.out
#SBATCH --output=slurm/TF_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

   IDN=${1}
 INPUT=${2}
OUTPUT=${3}

echo "TF: IDN:${IDN} + INPUT:${INPUT} -> OUTPUT:${OUTPUT}"

if [ "${INPUT}" == "" ] || [ ! -e ${INPUT} ]
then
	echo "TF: \"${INPUT}\" doesn't exist!"
#	scriptFailed "TF"
	exit 1
fi

if [ "${OUTPUT}" == "" ]
then
	echo "TF: Output not defined. Using input:${INPUT}"
	OUTPUT=${INPUT}
fi

LABEL=$(echo $(basename ${OUTPUT}) | tr '.' '_')

if [ -e ${LABEL}.transfer.done ]
then
	echo "TF: $LABEL already started!"
	exit 0
fi

if [ "${OUTPUT}" == "$(basename ${OUTPUT})" ]
then
	echo "TF: Output file has not pathing. Storing in temp folder."
	OUTPUT="~/projects/NeSI_Transfer/${IDN}/${OUTPUT}"
fi

echo "TF: ${INPUT} -> ${OUTPUT}"
date

CMD="ssh globus transfer --encrypt --perf-cc 4 --perf-p 8 --label \"${LABEL}\" -- nz#uoa/$(pwd)/${INPUT} nesi#otago-dtn01/${OUTPUT}"
echo "TF: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "TF:  ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "TF: ${OUTPUT} failed!"
#	scriptFailed "TF"
	exit 1
fi

touch ${LABEL}.transfer.done

#storeMetrics	# Unless we're going to wait around for this, why do we need the metrics?
