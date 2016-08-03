#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=2048
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/DC_%A_%a.out
#SBATCH --output=slurm/DC_%A_%a.out

source /projects/uoo00032/Resources/bin/baserefs.sh

PLATFORM=${1}
CONTIG=$([ "${2}" == "" ] && echo -ne "${CONTIGA[$SLURM_ARRAY_TASK_ID]}" || echo -ne "${2}")	# Input value or specific value.

   INPUT=printreads/${CONTIG}.bam
  OUTPUT=depth/${CONTIG}

echo "DC: ${INPUT} + ${PLATFORM} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "DC: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "DC"
	exit 1
fi

platformBED=${PLATFORMS}/${PLATFORM}.bed
  genderBED=${PLATFORMS}/$([ "$PLATFORM" == "Genomic" ] && echo -ne "AV5" || echo -ne "$PLATFORM" ).bed
  
# Special cases for X and Y depth of covereage as the X/YPAR1 and X/YPAR2 regions are distorted.
# Genomic Y is rife with repeat sequences that inflate coverage so use AV5 region for that.
# X: █▄▄▄▄▄▄▄█
# Y: _▄▄▄▄▄▄▄_
if [ "${CONTIG}" == "X" ]; then
	platformFile=${genderBED}
	actualContig=${TRUEX}
elif [ "${CONTIG}" == "Y" ]; then
	platformFile=${genderBED}
	actualContig=${TRUEY}
else
	platformFile=${platformBED}
	actualContig=${CONTIG}
fi


GATK_PROC=DepthOfCoverage
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
-L ${platformFile} \
-L ${actualContig} \
-isr INTERSECTION \
--omitLocusTable \
--omitDepthOutputAtEachBase \
--omitIntervalStatistics \
-nt ${SLURM_JOB_CPUS_PER_NODE}"


module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar $GATK ${GATK_ARGS} -I ${INPUT} -o ${OUTPUT}"
echo "DC: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "DC: ${INPUT} + ${PLATFORM} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "DC: ${OUTPUT} failed!"
#	scriptFailed "DC"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics
