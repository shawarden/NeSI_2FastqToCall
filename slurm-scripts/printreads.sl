#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=4092
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/PR_%j.out
#SBATCH --output=slurm/PR_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

CONTIG=${CONTIGA[$SLURM_ARRAY_TASK_ID]}
 INPUT=markdup/${CONTIG}.bam
  BQSR=baserecal/${CONTIG}.firstpass
OUTPUT=printreads/${CONTIG}.bam

echo "PR: ${INPUT} + ${BQSR} -> ${OUTPUT}"
date

if [ ! -e ${INPUT} ]; then
	echo "PR: Input file \"${INPUT}\" doesn't exist!"
#	scriptFailed "PR"
	exit 1
fi

GATK_PROC=PrintReads
GATK_ARGS="-T ${GATK_PROC} \
-R ${REFA} \
--bam_compression 9 \
-L ${CONTIG} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

module load ${MOD_JAVA}

CMD="$(which srun) $(which java) ${JAVA_ARGS} -jar ${GATK} ${GATK_ARGS} -I ${INPUT} -BQSR ${BQSR} -o ${OUTPUT}"
echo "PR: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "PR: ${INPUT} to ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "PR: ${INPUT} failed!"
#	scriptFailed "PR"
	exit 1
fi

rm ${INPUT} ${INPUT%.bam}.bai ${INPUT%.bam}.metrics

touch ${OUTPUT}.done

storeMetrics
