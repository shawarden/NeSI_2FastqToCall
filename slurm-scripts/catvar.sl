#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=16384
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/CV_%j.out
#SBATCH --output=slurm/CV_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

 INPUT=${1}
OUTPUT=${2}

echo "CV: ${INPUT} -> ${OUTPUT}"
date

mergeList=""
for file in ${INPUT}; do
	if [ -e $file ]; then
		mergeList="${mergeList} -V ${file}"
	else
		echo "CV: $file doesn't exist!"
#		scriptFailed "CarVar"
		exit 1
	fi
done

GATK_PROC=org.broadinstitute.gatk.tools.CatVariants
GATK_ARGS="${GATK_PROC} \
-R ${REFA} \
--assumeSorted"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -cp $GATK ${GATK_ARGS} ${mergeList} -out ${OUTPUT}"
echo "CV: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "CV: ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "CV: ${OUTPUT} failed!"
#	scriptFailed "CarVar"
	exit 1
fi

rm ${INPUT}

touch ${OUTPUT}.done

storeMetrics
