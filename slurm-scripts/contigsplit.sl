#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=2048
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/CS_%A_%a.out
#SBATCH --output=slurm/CS_%A_%a.out

source /projects/uoo00032/Resources/bin/baserefs.sh

SAMPLE=${1}
INPUT=${2}
BLOCK=${3}

CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
OUTPUT=split/${CONTIG}/split_${BLOCK}.bam


echo "CS: ${INPUT} by contig ${CONTIG} to ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "CS: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "CS"
	exit 1
fi

CMD="${SAMTOOLS} view -bh -@ 1 ${INPUT} ${CONTIG} > ${OUTPUT}"
echo "CS: ${CMD}" | tee -a ../commands.txt

eval ${CMD}
passed=$?

echo "CS: ${INPUT} to ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]
then
	echo "CS: ${INPUT} to ${CONTIG} Failed!"
#	scriptFailed "CS"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics run
