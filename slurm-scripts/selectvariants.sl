#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time 06:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/SV_%j.out
#SBATCH --output=slurm/SV_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

IDN=${1}
INPUT=${2}
OUTPUT=${INPUT%.bam}.vcf.gz

echo "FP: ${INPUT} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "FP: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "HC"
	exit 1
fi

GATK_PROC=SelectVariants
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-L ${COMMON} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar $GATK ${GATK_ARGS} -I ${INPUT} -o ${OUTPUT}"
echo "FP: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "FP: ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "FP: ${OUTPUT} failed!"
#	scriptFailed "HC"
	exit 1
fi

touch ${OUTPUT}.done

sbatch -J SV_${IDN} ${SLSBIN}/selectvariants.sl ${OUTPUT} ${3}

storeMetrics
