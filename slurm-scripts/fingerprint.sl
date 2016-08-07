#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time 06:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/FPH_%j.out
#SBATCH --output=slurm/FPH_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

IDN=${1}
INPUT=${2}
OUTPUT=${INPUT%.bam}.vcf.gz

HEADER="FHC"

echo "$HEADER: ${INPUT} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile; then exit 1; fi
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

GATK_PROC=HaplotypeCaller
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-L ${COMMON} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar $GATK ${GATK_ARGS} -I ${INPUT} -o ${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

${CMD}
if cmdFailed; then exit 1; fi

# Move output to final location
if ! finalOut; then exit 1; fi

touch ${OUTPUT}.done

sbatch -J FPS_${IDN} ${SLSBIN}/selectvariants.sl ${OUTPUT} ${3}

storeMetrics
