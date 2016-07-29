#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=02:00:00
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

CMD1="$(which srun) ${SAMTOOLS} cat -h ${HEADER} -o ${OUTPUT} ${INPUT}"
CMD2="$(which srun) ${SAMTOOLS} index ${OUTPUT} ${INDEX_OUTPUT}"
echo "CR: ${CMD1}, ${CMD2}" | tee -a commands.txt

${CMD1}
passed=$?

if [ $passed -eq 0 ]; then
	${CMD2}
	passed=$?
fi

echo "CR: ${INPUT} -> ${OUTPUT} & ${INDEX_OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "CR: ${INPUT} -> ${OUTPUT} & ${INDEX_OUTPUT} Failed!"
#	scriptFailed "CR"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics


