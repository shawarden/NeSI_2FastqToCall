#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=16384
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/CV_%j.out
#SBATCH --output=slurm/CV_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

 FILES=${1}
OUTPUT=${2}
IDN=$(echo $SLURM_JOB_NAME | cut -d'_' -f2)

HEADER="CV"

echo $HEADER $INPUT "->" $OUTPUT
date

mergeList=""
for INPUT in ${FILES}; do
	if ! inFile; then
		exit 1
	else 
		mergeList="${mergeList} -V ${INPUT}"
	fi
done

# Make sure input and target folders exists and that output file does not!
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

GATK_PROC=org.broadinstitute.gatk.tools.CatVariants
GATK_ARGS="${GATK_PROC} \
-R ${REFA} \
--assumeSorted"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -cp $GATK ${GATK_ARGS} ${mergeList} -out ${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 1
fi

# Move output to final location
if ! finalOut; then exit 1; fi

rm $FILES && echo "$HEADER: Purged input files!"

touch ${OUTPUT}.done

# Start transfers for variants file and index.
sbatch -J TV_${IDN} ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}
sbatch -J TVI_${IDN} ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}.tbi

storeMetrics
