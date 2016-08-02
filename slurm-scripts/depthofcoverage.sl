#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/DC_%j.out
#SBATCH --output=slurm/DC_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

  CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
   INPUT=printreads/${CONTIG}/printreads.bam
PLATFORM=${1}
  OUTPUT=depth/${CONTIG}/depth

echo "DC: ${INPUT} + ${PLATFORM} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "DC: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "DC"
	exit 1
fi

platformBED=${PLATFORMS}/${PLATFORM}.bed
  genderBED=${PLATFORMS}/$([ "${PLATFORM}" == "Genomic" ] && echo -ne "AV5" || echo -ne "${PLATFORM}" ).bed
  
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
-L ${PLATFORM} \
-L ${CONTIG} \
-isr INTERSECTION \
--omitDepthOutputAtEachBase \
--omitLocusTable \
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
