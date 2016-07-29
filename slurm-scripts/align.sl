#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/PA_%j.out
#SBATCH --output=slurm/PA_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

   SAMPLE=${1}
READGROUP=${2}
    BLOCK=${3}
   OUTPUT=${4}

READ1=blocks/${READGROUP}_R1_${BLOCK}.fastq.gz
READ2=blocks/${READGROUP}_R2_${BLOCK}.fastq.gz

echo "PA: ${READGROUP} ${READ1} ${READ2} to ${OUTPUT}"
date

if [ ! -e ${READ1} ] || [ ! -e ${READ2} ]; then
	echo "PA: A read files doesn't exist!"
#	scriptFailed "PA"
	exit 1
fi

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

echo "PA: ${BWA} mem -M -t ${SLURM_JOB_CPUS_PER_NODE} -R @RG'\t'$RG_ID'\t'$RG_PL'\t'$RG_PU'\t'$RG_LB'\t'$RG_SM $REF $READ1 $READ2 | ${SAMTOOLS} view -bh - > ${OUTPUT}" | tee -a ../commands.txt

${BWA} mem -M -t ${SLURM_JOB_CPUS_PER_NODE} -R @RG'\t'$RG_ID'\t'$RG_PL'\t'$RG_PU'\t'$RG_LB'\t'$RG_SM $REF $READ1 $READ2 | ${SAMTOOLS} view -bh - > ${OUTPUT}

passed=$((${PIPESTATUS[0]} + ${PIPESTATUS[1]}))

#mv ${JOB_TEMP_DIR}/${OUTPUT} ${OUTPUT} && rm ${JOB_TEMP_DIR}/${OUTPUT}

echo "PA: ${READGROUP} ${READ1} and ${READ2} ran for $(($SECONDS / 3600))h $((($SECONDS %3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ] || [ $(stat --printf="%s" ${OUTPUT}) -lt 50 ]; then
	echo "PA: ${READ1} and ${READ2} failed!"
#	scriptFailed "PA"
	exit 1
fi

rm ${READ1} ${READ2}

touch ${OUTPUT}.done

storeMetrics run
