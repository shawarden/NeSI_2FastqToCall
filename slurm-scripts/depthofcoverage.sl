#!/bin/bash
#SBATCH --job-name		DepthOfCoverage
#SBATCH --time			0-00:30:00
#SBATCH --mem-per-cpu	2048
#SBATCH --cpus-per-task	4
#SBATCH --constraint	avx
#SBATCH --array			1-84
#SBATCH --error			slurm/DC_%A_%a.out
#SBATCH --output		slurm/DC_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

PLATFORM=${1}
CONTIG=$([ "${2}" == "" ] && echo -ne "${CONTIGA[$SLURM_ARRAY_TASK_ID]}" || echo -ne "${2}")	# Input value or specific value.

 INPUT=printreads/${CONTIG}.bam
OUTPUT=depth/${CONTIG}

HEADER="DC"

echo "$HEADER: ${INPUT} + ${PLATFORM} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit 1; fi
if ! outFile; then exit 1; fi

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
echo "$HEADER: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics
