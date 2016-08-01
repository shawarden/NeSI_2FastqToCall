#!/bin/bash

#############################################################
# Generates job sequence for a given individual/sample list #
#############################################################

# Get base values
source /projects/uoo00032/Resources/bin/baserefs.sh

SAMPLE=${1}
READ1=${FASTQS}/${2}
READ2=${FASTQS}/${3}
PLATFORM=${4}
LOCATION=${5}

printf "%-12s%s\n" "SampleID" "${SAMPLE}"
printf "%-12s%s\n" "Read 1" "${READ1}"
printf "%-12s%s\n" "Read 2" "${READ2}"
printf "%-12s%s\n" "PLATFORM" "${PLATFORM}"

IDN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')
DNA=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $2}')
LIB=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $3}')
RUN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $4}')

# Platform setup.
platformBED=${PLATFORMS}/${PLATFORM}.bed
  genderBED=${PLATFORMS}/$([ "${PLATFORM}" == "Genomic" ] && echo "AV5" || echo "${PLATFORM}" ).bed
  genderSRC=${genderBED%.bed}.sh
  
   CUR_PATH=$(pwd)
if [ "$LOCATION" == "" ] || [ "$LOCATION" == "scratch" ]; then
	WORK_PATH=/scratch/jobs/$USER
else
	WORK_PATH=/projects/uoo00032/Alignments
fi


SAMPLE_PATH=${WORK_PATH}/${IDN}/${DNA}_${LIB}_${RUN}

mkdir -p ${SAMPLE_PATH}
mkdir -p ${SAMPLE_PATH}/slurm

date '+%Y%m%d_%H%M%S' > ${WORK_PATH}/${IDN}/starttime.txt

# Purge existing merge dependencies, just in case!
if [ -d ${SAMPLE_PATH}/mergeDeps ]; then
	printf "%-12s%s\n" "Purging" "mergeDep!"
	rm -r ${SAMPLE_PATH}/mergeDeps
fi

cd ${SAMPLE_PATH}

# Split read 1 and 2 into chunks
if [ ! -e ${SAMPLE}_R1_split.done ] || [ ! -e ${SAMPLE}_R2_split.done ]; then
	# One of the reads isn't finished yet. Run them both.
	
	if [ ! -e ${SAMPLE}_R1_split.done ]; then
		# Read 1 split isn't complete, run it now.
		read1Size=$(ls -lah ${READ1} | awk '{print $5}')
		DEP_SR1=$(sbatch -J RS_${SAMPLE}_R1_${read1Size} ${SLSBIN}/readsplit.sl ${SAMPLE} R1 ${READ1} | awk '{print $4}')
		printf "%-12s%s\n" "Split 1" "${DEP_SR1}"
	else
		printf "%-12s%s\n" "Split 1" "Done!"
	fi
	
	if [ ! -e ${SAMPLE}_R2_split.done ]; then
		# Read 2 split isn't complete. run it now.
		read2Size=$(ls -lah ${READ2} | awk '{print $5}')
		DEP_SR2=$(sbatch -J RS_${SAMPLE}_R2_${read2Size} ${SLSBIN}/readsplit.sl ${SAMPLE} R2 ${READ2} | awk '{print $4}')
		printf "%-12s%s\n" "Split 2" "${DEP_SR2}"
	else
		printf "%-12s%s\n" "Split 2" "Done!"
	fi
else
	printf "%-12s%s\n" "Split" "Done! Aligning..."
	# Both reads are done. Build read list and spool alignments.
	READGROUP=$(cat blocks/R1_ReadGroup.txt)
	readBlocks=$(($(find ./blocks -type f -iname "${READGROUP}_R1_*.fastq.gz.done" | wc -l) -1))
	for i in $(seq 0 ${readBlocks}); do
		if [ $i -lt $readBlocks ]; then
			nextBlock=$(($i + 1))
		else
			nextBlock=$i
		fi
		${PBIN}/check_blocks.sh ${SAMPLE} ${READGROUP} R1 $(printf '%05d' $i) $(printf '%05d' $nextBlock)
	done
fi

exit 0
