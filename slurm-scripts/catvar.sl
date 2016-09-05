#!/bin/bash
#SBATCH --job-name		CatVariants
#SBATCH --time			0-03:00:00
#SBATCH --mem-per-cpu	16384
#SBATCH --cpus-per-task	1
#SBATCH --constraint	avx
#SBATCH --error			slurm/CV_%j.out
#SBATCH --output		slurm/CV_%j.out

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
		exit $EXIT_IO
	else 
		mergeList="${mergeList} -V ${INPUT}"
	fi
done

# Make sure input and target folders exists and that output file does not!
if ! outDirs; then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

GATK_PROC=org.broadinstitute.gatk.tools.CatVariants
GATK_ARGS="${GATK_PROC} \
-R ${REFA} \
--assumeSorted"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -cp $GATK ${GATK_ARGS} ${mergeList} -out ${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit $EXIT_PR
fi

# Move output to final location
if ! finalOut; then exit $EXIT_MV; fi

rm $FILES && echo "$HEADER: Purged input files!"

touch ${OUTPUT}.done

# Start transfers for variants file and index.
if ! ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}; then
	echo "$HEADER: Transfer failed!"
	exit $EXIT_TF
fi

if ! ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}.tbi; then
	echo "$HEADER: Transfer index failed!"
	exit $EXIT_TF
fi
