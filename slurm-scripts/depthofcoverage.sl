#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/depthofcoverage_%j.out
#SBATCH --output=slurm/depthofcoverage_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

   INPUT=${1}
  CONTIG=${2}
PLATFORM=${3}
  OUTPUT=${4}

echo "Depth: ${INPUT} + ${PLATFORM} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "Depth: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "Depth"
	exit 1
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
echo "Depth: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "Depth: ${INPUT} + ${PLATFORM} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "Depth: ${OUTPUT} failed!"
#	scriptFailed "Depth"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics
