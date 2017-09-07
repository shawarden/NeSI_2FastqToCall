#!/bin/bash
#SBATCH --job-name		CatReads
#SBATCH --time			60
#SBATCH --mem-per-cpu	4096
#SBATCH --cpus-per-task	1
#SBATCH --constraint	avx
#SBATCH --error			slurm/CR_%j.out
#SBATCH --output		slurm/CR_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

FILES=${1}
OUTPUT=${2}
IDN=$(echo $SLURM_JOB_NAME | cut -d'_' -f2)

BAMHEAD=$(echo ${FILES} | awk '{print $1}')

HEADER="CR"

# Get local node freespace for /tmp and /scratch
DF_OUT=$(df -a)
DF_TMP=$(echo "$DF_OUT" | grep /tmp)
DF_SCRATCH=$(echo "$DF_OUT" | grep /scratch)

# If node's temp folder has enough space, write output locally, otherwise leave at main destination.
if $(echo $DF_TMP | awk '{if (($4/1024/1024) > 222.0) {exit 0} else {exit 1}}'); then
	echo "$HEADER: Writing to local node /tmp folder. $DF_TMP"
	OUTDIR=$JOB_TEMP_DIR
elif $(echo $DF_SCRATCH | awk '{if (($4/1024/1024) > 222.0) {exit 0} else {exit 1}}'); then
	echo "$HEADER: Not enough space on local node. Writing to scratch. $DF_SCRATCH"
	OUTDIR=$SCRATCH_DIR
else
	echo "$HEADER: No enough space on local node or scratch disk for write. Writing to final destination."
	df -ah
fi

echo $HEADER: $FILES + "Header($BAMHEAD) ->" $OUTPUT
date

for INPUT in ${FILES}; do
	if ! inFile; then exit $EXIT_IO; fi
done

if ! outDirs; then exit $EXIT_IO; fi
if ! outFile; then exit $EXIT_IO; fi

CMD="srun ${SAMTOOLS} cat -h ${BAMHEAD} -o ${OUTDIR}/${OUTPUT} ${FILES}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

JOBSTEP=0

if ! ${CMD}; then
	cmdFailed $?
	exit ${JOBSTEP}${EXIT_PR}
fi

storeMetrics

# Move output to final location
if [ "$OUTDIR" == "$JOB_TEMP_DIR" ]; then
	if ! finalOut; then exit $EXIT_MV; fi
elif [ "$OUTDIR" == "$SCRATCH_DIR" ]; then
	if ! scratchOut; then exit $EXIT_MV; fi
fi

touch ${OUTPUT}.done

sbatch $(dispatch "RI") -J RI_${IDN} ${SLSBIN}/catreadsindex.sl ${OUTPUT} ${OUTPUT%.bam}.bai

if ! . ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}; then
	echo "$HEADER: Transfer failed!"
	exit $EXIT_TF
fi
