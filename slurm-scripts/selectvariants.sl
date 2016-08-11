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
if ! inFile;  then exit 1; fi
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

GATK_PROC=SelectVariants
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-L ${COMMON} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar $GATK ${GATK_ARGS} -I ${INPUT} -o ${JOB_TEMP_DIR}/${OUTPUT}"
echo "SV: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 1
fi

# Move output to final location
if ! finalOut; then exit 1; fi

touch ${OUTPUT}.done

sbatch -J SFP_${IDN} ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}

storeMetrics
