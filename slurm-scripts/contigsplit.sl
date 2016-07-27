#!/bin/bash
#SBATCH -A uoo00032
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=2048
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/contigsplit_%j.out
#SBATCH --output=slurm/contigsplit_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

SAMPLE=${1}
CONTIG=${2}
 INPUT=${3}
OUTPUT=${4}

echo "ContigSplit: ${INPUT} by contig ${CONTIG} to ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "ContigSplit: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "ContigSplit"
	exit 1
fi

CMD="${SAMTOOLS} view -bh -@ 1 ${INPUT} ${CONTIG}"
echo "ContigSplit: ${CMD} > ${OUTPUT}" | tee -a ../commands.txt

${CMD} > ${OUTPUT}
passed=$?

echo "ContigSplit: ${INPUT} to ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]
then
	echo "ContigSplit: ${INPUT} to ${CONTIG} Failed!"
#	scriptFailed "ContigSplit"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics run
