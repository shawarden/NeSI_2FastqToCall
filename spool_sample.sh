#!/bin/bash

#############################################################
# Generates job sequence for a given individual/sample list #
#############################################################

# Get base values
source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

SAMPLE=${1}
READ1=${FASTQS}/${2}
READ2=${FASTQS}/${3}
PLATFORM=${4}
LOCATION=${5}

printf "%-9s%s\n" "SampleID" "${SAMPLE}"
printf "%-9s%s\n" "Read 1" "${READ1}"
printf "%-9s%s\n" "Read 2" "${READ2}"
printf "%-9s%s\n" "PLATFORM" "${PLATFORM}"

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
	printf "%-9s%s\n" "Purging" "mergeDep!"
	rm -r ${SAMPLE_PATH}/mergeDeps
fi

cd ${SAMPLE_PATH}


# Dispatch alignemnt & sort arrays.
# Alignemnt array doens't have an dependency yet since ReadSplit needs to know what the aligner's JobID is to update it.
printf "%-9s " "Dispatch (Align"
DEP_PA=$(sbatch -J PA_${SAMPLE} --begin=now+10hour --array=0-999 ${SLSBIN}/alignar.sl ${SAMPLE} ${READGROUP} | awk '{print $4}')
if [ $? -ne 0 ] || [ "$DEP_PA" == "" ]; then
	printf "FAILED!"
	exit 1
else
	printf "%s " "${DEP_PA}"
fi

# ReadSplit also needs to know Sort's JobID for the same reason.
printf "%s " "-> Sorter"
DEP_SS=$(sbatch -J SS_${SAMPLE} $(depCheck ${DEP_PA}) --array=0-999 ${SLSBIN}/sortar.sl ${SAMPLE} ${READGROUP} | awk '{print $4}')
if [ $? -ne 0 ] || [ "$DEP_SS" == "" ]; then
	printf "FAILED!"
	exit 1
else
	printf "%s) " "${DEP_SS}"
fi

printf "%s " "<- Split Reads"

# Split read 1 and 2 into chunks
splitReadArray=""

if [ ! -e ${SAMPLE}_R1_split.done ]; then
	splitReadArray=$(appendList "$splitReadArray" 1 ",")
fi

if [ ! -e ${SAMPLE}_R2_split.done ]; then
	splitReadArray=$(appendList "$splitReadArray" 2 ",")
fi

# Read 1 split isn't complete, run it now.
read1Size=$(ls -lah ${READ1} | awk '{print $5}')
read2Size=$(ls -lah ${READ2} | awk '{print $5}')

if [ "$splitReadArray" != "" ]; then
	DEP_SR=$(sbatch -J RS_${SAMPLE}_${read1Size}_${read2Size} -a $splitReadArray ${SLSBIN}/readsplit.sl ${SAMPLE} ${READ1} ${READ2} ${DEP_PA} ${DEP_SS} | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_SR" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "${DEP_SR}"
		
		# Now that we have the ReadSplit jobID, assign the entire Alignment array to it.
		scontrol update JobID=${DEP_PA} StartTime=now Dependency=afterok:${DEP_SR}
	fi
else
	printf "done -> Aligning..."
	# Both reads are done. Build read list and spool alignments.
	READGROUP=$(cat blocks/R1_ReadGroup.txt)
	readBlocks=$(($(find ./blocks -type f -iname "${READGROUP}_R1_*.fastq.gz.done" | wc -l) -1))
	for i in $(seq 0 ${readBlocks}); do
		if [ $i -lt $readBlocks ]; then
			nextBlock=$(($i + 1))
		else
			nextBlock=$i
		fi
		${PBIN}/check_blocks.sh ${SAMPLE} ${READGROUP} R1 $(printf '%05d' $i) $(printf '%05d' $nextBlock $DEP_PA $DEP_SS)
	done
fi




