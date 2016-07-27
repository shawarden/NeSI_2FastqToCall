#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/mergecontigs_%j.out
#SBATCH --output=slurm/mergecontigs_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

CONTIG=${1}
OUTPUT=${2}
 CLEAN=${3}

echo "MergeSAM: ${CONTIG} -> ${OUTPUT}"
date

contigMerBlocks=$(find . -type f -iwholename "*/contig/${CONTIG}/*.bam" -printf '%h\0%d\0%p\n' | sort -t '\0' -n | awk -F '\0' '{print $3}')
if [ "${contigMerBlocks}" == "" ]; then
	echo "MergeSAM: Merge contig ${CONTIG} contains no files!"
#	scriptFailed "MergeSAM"
	exit 1
else
	echo "MergeSAM: Merge contig ${CONTIG} will run ${contigMerBlocks}"
fi

mergeList=""
for file in ${contigMerBlocks}; do
	if [ -e $file ]; then
		mergeList="${mergeList} I=${file}"
	else
		echo "MergeSAM: Input file $file doesn't exist!"
#		scriptFailed "MergeSAM"
		exit 1
	fi
done

if [ "$mergeList" == "" ]; then
	echo "MergeSAM: No inputs defined!"
#	scriptFailed "MergeSAM"
	exit 1
fi

PIC_ARGS="SORT_ORDER=coordinate \
USE_THREADING=true \
CREATE_INDEX=true \
COMPRESSION_LEVEL=9 \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
TMP_DIR=${JOB_TEMP_DIR}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${PICARD} MergeSamFiles ${PIC_ARGS} ${mergeList} O=${OUTPUT}"
echo "MergeSAM: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

# Report run time.
echo "MergeSAM: ran for $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "MergeSAM: ${CONTIG} failed!"
#	scriptFailed "MergeSAM"
	exit 1
fi

rm ${contigMerBlocks}

touch ${OUTPUT}.done

storeMetrics
