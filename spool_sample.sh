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
MULTI_RUN=${6}

if [ "$LOCATION" == "scratch" ]; then
	WORK_PATH="/scratch/jobs/$USER"
else
	WORK_PATH="/projects/uoo00032/Alignments"
fi

printf "%-22s%s\n" "SampleID" "${SAMPLE}"
printf "%-22s%s\n" "Read 1" "${READ1}"
printf "%-22s%s\n" "Read 2" "${READ2}"
printf "%-22s%s\n" "Platform" "${PLATFORM}"
printf "%-22s%s\n" "Location" "${WORK_PATH}"
[ "$MULTI_RUN" == "" ] && printf "%-22s%s\n" "MutliRun" "YES"

IDN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')
DNA=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $2}')
LIB=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $3}')
RUN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $4}')

INDIVIDUAL_PATH=${WORK_PATH}/${IDN}
SAMPLE_PATH=${INDIVIDUAL_PATH}/${DNA}_${LIB}_${RUN}

if ! mkdir -p ${SAMPLE_PATH}/slurm; then
	echo "Error creating output folder!"
	exit 1
fi

if ! mkdir -p ${INDIVIDUAL_PATH}/slurm; then
	echo "Error creating output folder!"
	exit 1
fi

printf "%-22s" "Command"
echo $0 ${@} | tee ${INDIVIDUAL_PATH}/jobReSubmit.sh
chmod +x ${INDIVIDUAL_PATH}/jobReSubmit.sh

date '+%Y%m%d_%H%M%S' >> ${WORK_PATH}/${IDN}/starttime.txt

cd ${SAMPLE_PATH}

printf "%-22s" "ReadSplitter"

##################################
# Split read 1 and 2 into chunks #
##################################

splitReadArray=""

for i in $(seq 1 2); do
	# Cycle through reads.
	if [ ! -e ${SAMPLE}_R${i}_split.done ]; then
		# Read# split isn't complete. Add to array.
		splitReadArray=$(appendList "$splitReadArray" ${i} ",")
	fi
done

############################
# Launch needed split jobs #
############################

if [ "$splitReadArray" != "" ]; then
	
	if [ ! -e $READ1 ] || [ ! -e $READ2 ]
	then
		echo "Read files may not exist!"
		exit 1
	fi
	
	####################
	# Get a fancy size #
	####################
	
	sizeString=" kMGTEPYZ"
	sizeBlock=0
	readSize=$(($(ls -la $READ1 | awk '{print $5}') + $(ls -la $READ2 | awk '{print $5}')))
	while [ $(echo "$readSize / 1024 > 0" | bc) -eq 1 ]; do
		#printf "%-12s %.0f%s\n" "Read size" $readSize $(echo ${sizeString:${sizeBlock}:1}Bytes | sed -e 's/ //g')
		readSize=$(echo "$readSize / 1024" | bc -l)
		sizeBlock=$((sizeBlock+1))
	done
	readSize=$(echo $(printf "%.0f" $readSize)${sizeString:${sizeBlock}:1}B | sed -e 's/ //g')
	
	# Split array contains data so run the missing split function.
	DEP_RS=$(sbatch $(dispatch "RS") -J RS_${SAMPLE}_${readSize} -a $splitReadArray $SLSBIN/readsplit.sl $SAMPLE $READ1 $READ2 $PLATFORM ${MULTI_RUN} | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_RS" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%sx%-1d [%s]\n" "${DEP_RS}" $(splitByChar "$splitReadArray" "," | wc -w) "$splitReadArray"
		echo $DEP_RS > ../lastJob.txt
	fi
else
	printf "done\n"
	if ! ${PBIN}/spool_merge.sh ${SAMPLE} ${PLATFORM} ${MULTI_RUN}; then
		cmdFailed $?
		exit $EXIT_PR
	fi
fi
