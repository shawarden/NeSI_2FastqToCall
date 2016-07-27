#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time 03:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/haplotypecaller_%j.out
#SBATCH --output=slurm/haplotypecaller_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

 INPUT=${1}
CONTIG=${2}
OUTPUT=${3}

echo "Haplo: Command-line: ${0} ${@}"

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

echo "Haplo: ${INPUT} + ${CONTIG} + ${intervalPloidy} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "Haplo: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "Haplo"
	exit 1
fi

GATK_PROC=HaplotypeCaller
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-L ${CONTIG} \
--sample_ploidy ${intervalPloidy} \
--emitRefConfidence GVCF \
--dbsnp ${DBSNP} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar $GATK ${GATK_ARGS} -I ${INPUT} -o ${OUTPUT}"
echo "Haplo: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "Haplo: ${INPUT} + ${CONTIG} + ${intervalPloidy} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "Haplo: ${OUTPUT} failed!"
#	scriptFailed "Haplo"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics
