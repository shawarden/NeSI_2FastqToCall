#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --time			0-03:00:00
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	8
#SBATCH --constraint	avx
#SBATCH --array			1-84
#SBATCH --error			slurm/HC_%A_%a.out
#SBATCH --output		slurm/HC_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

CONTIG=$([ "$1" != "" ] && echo -ne "$1" || echo -ne "${CONTIGARRAY[$SLURM_ARRAY_TASK_ID]}")

HEADER="HC"

 INPUT=printreads/${CONTIG%:*}.bam	# Strip any contig coordinates.
OUTPUT=haplo/${CONTIG}.g.vcf.gz

if [ -e coverage.sh ]; then
	# Gender file exists so obtain values from it.
	source coverage.sh
fi

# Calculate ploidy based on number of X and Y chromosomes.
if [ "$CONTIG" == "${XPAR1}" ] || [ "$CONTIG" == "${XPAR2}" ]; then
	# XPAR1 and XPAR2 ploidies are the number of X and Y chromosomes combines. 1+1=2, 2+0=2, etc.
	intervalPloidy=$(($X_CHROMOSOMES + $Y_CHROMOSOMES))
elif [ "$CONTIG" == "${TRUEX}" ]; then
	# TRUEX ploidy is the number of X chromosomes. 1, 2, etc. 0 should never happen as the gender determination job should fail on that.
	intervalPloidy=$X_CHROMOSOMES
elif [ "$CONTIG" == "Y" ] || [ "$CONTIG" == "${TRUEY}" ]; then
	# Y or TRUEY ploidy is the number of Y chromosomes or 1, whichever is higher. so 0=1, 1=1, 2=2, etc.
	[ $Y_CHROMOSOMES -gt 1 ] && intervalPloidy=$Y_CHROMOSOMES || intervalPloidy=1
else
	# Non Gender chromosome ploidy is 2.
	intervalPloidy=2
fi

echo "$HEADER: ${INPUT} + ${CONTIG}c + ${intervalPloidy}p -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit 10; fi
if ! outDirs; then exit 10; fi
if ! outFile; then exit 10; fi

GATK_PROC=HaplotypeCaller
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-L ${CONTIG} \
--sample_ploidy ${intervalPloidy} \
--emitRefConfidence GVCF \
--dbsnp ${DBSNP} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar $GATK ${GATK_ARGS} -I ${INPUT} -o ${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 15
fi

# Move output to final location
if ! finalOut; then exit 20; fi

touch ${OUTPUT}.done

storeMetrics
