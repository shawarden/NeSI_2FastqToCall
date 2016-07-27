#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/catreads_%j.out
#SBATCH --output=slurm/catreads_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

INPUT=${1}
OUTPUT=${2}

HEADER=$(echo ${INPUT} | awk '{print $1}')

echo "CatReads: ${INPUT} + HeaderOf(${HEADER})-> ${OUTPUT}"
date

for file in ${INPUT}; do
	if [ ! -e ${file} ]; then
		echo "CatReads: \"${file}\" doesn't exist!"
#		scriptFailed "CarReads"
		exit 1
	fi
done

inputCount=$(echo ${INPUT} | wc -w)
contigCount=$(echo ${CONTIGS} | wc -w)

if [ $inputCount -ne $contigCount ]; then
	echo "CatReads: Input file count (${inputCount}) doesn't match contig count (${contigCount})!"
#	scriptFailed "CarReads"
	exit 1
else
	echo "CatReads: Input file count (${inputCount}) matches contig count (${contigCount})"
fi

CMD="$(which srun) ${SAMTOOLS} cat -h ${HEADER} -o ${OUTPUT} ${INPUT}"
echo "CatReads: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "CatReads: ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "CatReads: ${INPUT} -> ${OUTPUT} Failed!"
#	scriptFailed "CatReads"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics


