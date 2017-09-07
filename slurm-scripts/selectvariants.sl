#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		SelectVariants
#SBATCH --time			0-06:00:00
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	8
#SBATCH --constraint	avx
#SBATCH --error			slurm/SV_%j.out
#SBATCH --output		slurm/SV_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

IDN=${1}
INPUT=${2}
OUTPUT=${IDN}_fingerprint.vcf.gz

echo "SV: ${INPUT} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit $EXIT_IO; fi
if ! outDirs; then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

GATK_PROC=SelectVariants
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-L ${COMMON} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="srun $(which java) ${JAVA_ARGS} -jar $GATK ${GATK_ARGS} -I ${INPUT} -o ${JOB_TEMP_DIR}/${OUTPUT}"
echo "SV: ${CMD}" | tee -a commands.txt

JOBSTEP=0

if ! ${CMD}; then
	cmdFailed $?
	exit ${JOBSTEP}${EXIT_PR}
fi

storeMetrics

# Move output to final location
if ! finalOut; then exit $EXIT_MV; fi

touch ${OUTPUT}.done

${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}

