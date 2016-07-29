#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=01:30:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/RI_%j.out
#SBATCH --output=slurm/RI_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

INPUT="${1}"
OUTPUT="${2}"

if [ "$OUTPUT" == "" ]; then
	echo "RI: No output defined. Usuing input."
	OUTPUT=${INPUT%.bam}.bai
fi

echo "RI: ${INPUT} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "RI: \"${file}\" doesn't exist!"
#	scriptFailed "CarReadsIndex"
	exit 1
fi

CMD="$(which srun) ${SAMTOOLS} index ${INPUT} ${OUTPUT}"
echo "RI: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "RI: ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "RI: ${INPUT} -> ${OUTPUT} failed!"
#	scriptFailed "RI"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics


