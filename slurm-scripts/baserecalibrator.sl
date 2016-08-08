#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/BR_%A_%a.out
#SBATCH --output=slurm/BR_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
 INPUT=markdup/${CONTIG}.bam
OUTPUT=baserecal/${CONTIG}.firstpass

HEADER="BR"

echo "$HEADER: ${INPUT} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit 1; fi
if ! outFile; then exit 1; fi

GATK_PROC=BaseRecalibrator
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-knownSites ${DBSNP} \
-knownSites ${MILLS} \
-L ${CONTIG} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${GATK} ${GATK_ARGS} -I ${INPUT} -o ${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics
