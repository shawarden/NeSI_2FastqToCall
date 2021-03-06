#!/bin/bash
#SBATCH --job-name		BlockAlign
#SBATCH --time			359
#SBATCH --mem-per-cpu	2048
#SBATCH --cpus-per-task	8
#SBATCH --constraint	avx
#SBATCH --array			0-999
#SBATCH --error			slurm/BA_%A_%a.out
#SBATCH --output		slurm/BA_%A_%a.out

export RAMDISK=2

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

SAMPLE=${1}

# Block depends on input or array id.
export BLOCK=$([ "$3" != "" ] && printf "%0${FASTQ_MAXZPAD}d" $3 || printf "%0${FASTQ_MAXZPAD}d" $SLURM_ARRAY_TASK_ID)
   
OUTPUT=split/${BLOCK}/contig_split
mkdir -p $(dirname ${OUTPUT})

export HEADER="BA"
   
READ1=blocks/R1_${BLOCK}.fastq.gz
READ2=blocks/R2_${BLOCK}.fastq.gz

READGROUP=$($CAT_CMD $READ1 | head -1 | awk -F'[@:]' '{print $2"_"$3"_"$4"_"$5"_"$10}' )

echo "$HEADER: $READGROUP $BLOCK $READ1 $READ2 -> $OUTPUT"
jobStats
date

if [ $(echo "$READGROUP" | wc -w) -gt 1 ]; then
	echo "$HEADER: Too many read-groups!"
	exit 1
fi

# Make sure input and target folders exists and that output file does not!
if ! (INPUT=${READ1}; inFile); then exit $EXIT_IO; fi
if ! (INPUT=${READ2}; inFile); then exit $EXIT_IO; fi
if ! outDirs; then exit $EXIT_IO; fi
#if ! outFile; then exit $EXIT_IO; fi

if [ ! -e sorted/${BLOCK}.done ]; then
	# Get readgroup blocks from either INFO_INFO_INFO_.. or INFO INFO INFO ...
	     INTRUMENT=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $1}')
	INSTRUMENT_RUN=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $2}')
	     FLOW_CELL=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $3}')
	     CELL_LANE=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $4}')
	         INDEX=$(echo ${READGROUP} | awk -F'[[:blank:]_]' '{print $5}')
	
	RG_ID="ID:${INTRUMENT}_${INSTRUMENT_RUN}_${FLOW_CELL}_${CELL_LANE}_${INDEX}"
	RG_PL="PL:Illumina"
	RG_PU="PU:${FLOW_CELL}.${CELL_LANE}"
	RG_LB="LB:${SAMPLE}"
	RG_SM="SM:$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')"
	
	HEADER="PA"
	JOBSTEP=0
	echo "$HEADER: Aligning!"
	
	echo $REFA | grep 38 && BWA_REF=$REFA || BWA_REF=$REF
	
	# Pipe output from alignment into sortsam
	CMD="srun ${BWA} mem -M -t ${SLURM_JOB_CPUS_PER_NODE} -R @RG'\t'$RG_ID'\t'$RG_PL'\t'$RG_PU'\t'$RG_LB'\t'$RG_SM $REF $READ1 $READ2 | ${SAMTOOLS} view -bh - > $SHM_DIR/align_${BLOCK}.bam"
	echo "$HEADER: ${CMD}" | tee -a ../commands.txt
	
	if ! eval ${CMD}; then
		cmdFailed $?
		exit ${JOBSTEP}${EXIT_PR}
	fi
	
	storeMetrics
	
	PA_SECONDS=$SECONDS
	SECONDS=0
	
	HEADER="SS"
	JOBSTEP=1
	
	module load ${MOD_JAVA}
	
	CMD="srun $(which java) ${JAVA_ARGS} -jar ${PICARD} SortSam ${PIC_ARGS} ${SORT_ARGS} INPUT=$SHM_DIR/align_${BLOCK}.bam OUTPUT=$SHM_DIR/sorted_${BLOCK}.bam"
	echo "$HEADER: ${CMD}" | tee -a ../commands.txt
	
	if ! ${CMD}; then
		cmdFailed $?
		exit ${JOBSTEP}${EXIT_PR}
	fi
	
	storeMetrics
	
	rm $SHM_DIR/align_${BLOCK}.bam && echo "$HEADER: Purged aligned block: $SHM_DIR/align_${BLOCK}.bam"
else
	echo "$HEADER: Alignment already completed!"
	SS_SECONDS=$SECONDS
	SECONDS=0
fi

HEADER="CS"
JOBSTEP=2
echo "$HEADER: Splitting by contig"
# Run multiple samview splits at once. capped at cpus assigned. Initially slow but should speed up with time.
if ! for contig in ${CONTIGARRAY[@]}; do echo $contig; done | (
	srun xargs -I{} --max-procs ${SLURM_JOB_CPUS_PER_NODE} bash -c '{
		OUTPUT="split/${BLOCK}/{}.bam"
		CMD="${SAMTOOLS} view -bh -o ${OUTPUT} $SHM_DIR/sorted_${BLOCK}.bam {}"
		echo "$HEADER: ${CMD}" | tee -a ../commands.txt
		${CMD}
		echo "result: $?"
	}'
); then
	cmdFailed $?
	exit ${JOBSTEP}${EXIT_PR}
fi

storeMetrics
JOBSTEP=""
SECONDS=$(($SECONDS + $SS_SECONDS + $PA_SECONDS))

# Put this after $CMD execution to remove empty files.
#		if ls ${OUTPUT} | awk "{print \$5}" | grep 37695; then
#			rm ${OUTPUT} && echo "Deleting ${OUTPUT} as no data for this contig."
#		fi

df -ah $SHM_DIR
rm $SHM_DIR/sorted_${BLOCK}.bam && echo "$HEADER: Purged sorted block: $SHM_DIR/sorted_${BLOCK}.bam"

# Remove input files.
#rm ${READ1} ${READ2} && echo "$HEADER: Purged read files!"

# Indicate completion.
touch ${OUTPUT}.done
