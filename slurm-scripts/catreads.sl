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

source /projects/uoo00032/Resources/bin/baserefs.sh

INPUT=${1}
OUTPUT=${2}
INDEX_OUTPUT=${OUTPUT%.bam}.bai
IDN=$(echo $SLURM_JOB_NAME | cut -d'_' -f2)

HEADER=$(echo ${INPUT} | awk '{print $1}')

echo "CR: ${INPUT} + HeaderOf(${HEADER})-> ${OUTPUT}"
date

for file in ${INPUT}; do
	if [ ! -e ${file} ]; then
		echo "CR: \"${file}\" doesn't exist!"
#		scriptFailed "CarReads"
		exit 1
	fi
done

inputCount=$(echo ${INPUT} | wc -w)
contigCount=$(echo ${CONTIGS} | wc -w)

if [ $inputCount -ne $contigCount ]; then
	echo "CR: Input file count (${inputCount}) doesn't match contig count (${contigCount})!"
#	scriptFailed "CarReads"
	exit 1
else
	echo "CR: Input file count (${inputCount}) matches contig count (${contigCount})"
fi

CMD="$(which srun) ${SAMTOOLS} cat -h ${HEADER} -o ${OUTPUT} ${INPUT}"
echo "CR: ${CMD1}, ${CMD2}" | tee -a commands.txt

${CMD}
passed=$?

echo "CR: ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "CR: ${INPUT} -> ${OUTPUT} Failed!"
#	scriptFailed "CR"
	exit 1
fi

touch ${OUTPUT}.done

sbatch -J TR_${IDN} ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}
sbatch -J RI_${IDN} ${SLSBIN}/catreadsindex.sl ${OUTPUT} ${OUTPUT%.bam}.bai

storeMetrics


