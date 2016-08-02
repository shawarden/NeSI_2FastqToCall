#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/BR_%j.out
#SBATCH --output=slurm/BR_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
 INPUT=markdup/${CONTIG}.bam
OUTPUT=baserecal/${CONTIG}.firstpass

echo "BR: ${INPUT} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "BR: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "BR"
	exit 1
fi

GATK_PROC=BaseRecalibrator
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-knownSites ${DBSNP} \
-knownSites ${MILLS} \
-L ${CONTIG} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${GATK} ${GATK_ARGS} -I ${INPUT} -o ${OUTPUT}"
echo "BR: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "BR: ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "BR: ${INPUT} Failed."
#	scriptFailed "BR"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics
