#!/bin/bash

#############################################################
# Generates job sequence for a given individual/sample list #
#############################################################

# Get base values
source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

function usage {
cat << EOF

*************************************
* This script spool up an alignment *
* run for the specified patient ID  *
* on the NeSI cluster.              *
*************************************
*
* usage: $0 options:
*
*********************************
*
* Required:
*   -i [IDLR]      Sample ID string: ID_DNA_LIB_RUN.
*   -r [FILE,FILE] Read files separated by comma.
*   -p [PLATFORM]  Capture platform/Exome chip.
*
* Options:
*   -e [step]      Entry point for script
*   -m             This sample has multiple runs
*                  and this is NOT the last run.
*                  Omit this option on final run.
*   -s             Run in scratch area?
*   -c [XY]        Chromosomal Gender: XX, XY, X0, XXX, etc.
*   -g [Gender]    Gender: Male/Female/Unknown
*
*********************************
EOF
}

#ssh nesipan "source .bash_profile; spool_sample.sh -i ${SAMPLE} -r $(basename ${READ1}),$(basename ${READ2}) -p ${PLATFORM} -s ${USE_SCRATCH} $( [ ! -z ${MULTI_RUN} ] && echo -ne "-m" ) $([ ! -z ${GENDER} ] && echo -ne "-g ${GENDER}") $([ ! -z $SEXCHR ] && echo -ne "-c $SEXCHR")"

while getopts "msc:g:i:p:r:" OPTION
do
	FILE=
	case $OPTION in
		m)
			export MULTI_RUN=="-m"
#			(echo "multirun enabled" 1>&2)
			;;
		s)
			export USE_SCRATCH=YES
#			(echo "using scratch" 1>&2)
			;;
		c)
			export SEXCHR=${OPTARG}
#			(echo "assuming sexchroms $SEXCHR" 1>&2)
			;;
		g)
			export GENDER=${OPTARG}
#			(echo "assuming gender $GENDER" 1>&2)
			;;
		i)
			export SAMPLE=${OPTARG}
#			(echo "sample $SAMPLE" 1>&2)
			;;
		p)
			export PLATFORM=${OPTARG}
#			(echo "platform $PLATFORM" 1>&2)
			;;
		r)
			oldIFS=$IFS
			IFS=','
			for file in ${OPTARG}; do
				if [ ! -e ${FASTQS}/$file ]; then
					echo "WARN: Reads file ${FASTQS}/$file does not exist."
					#exit 1
#				else
#					(echo "read: ${FASTQS}/${file}" 1>&2)
				fi
			done
			READ1=${FASTQS}/$(echo ${OPTARG} | awk '{print $1}')
			READ2=${FASTQS}/$(echo ${OPTARG} | awk '{print $2}')
			IFS=$oldIFS
			;;
		?)
			echo "FAILURE: ${OPTION} ${OPTARG} is not valid!"
			usage
			exit 1
			;;
	esac
done

if [ "${READ1}" == "" ] || [ "${READ2}" == "" ] || [ "${SAMPLE}" == "" ] || [ "${PLATFORM}" == "" ]; then
	usage
	exit 1
fi

if [ "$USE_SCRATCH" == "YES" ]; then
	WORK_PATH="/scratch/jobs/$USER"
else
	WORK_PATH="/projects/uoo00032/Alignments"
fi

printf "%-22s%s\n" "SampleID" "${SAMPLE}"
printf "%-22s%s\n" "Read 1" "${READ1}"
printf "%-22s%s\n" "Read 2" "${READ2}"
printf "%-22s%s\n" "Platform" "${PLATFORM}"
printf "%-22s%s\n" "Location" "${WORK_PATH}"
[ "$MULTI_RUN" != "" ] && printf "%-22s%s\n" "MutliRun" "YES"
[ "$GENDER" != "" ] && printf "%-22s%s\n" "Gender" "$GENDER"
[ "$SEXCHR" != "" ] && printf "%-22s%s\n" "Sex Chroms" "$SEXCHR"

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

cd ${INDIVIDUAL_PATH}

[ "$MULTI_RUN" != "" ] && echo "${SAMPLE}" > multirun.txt
[ "$SEXCHR" != "" ] && echo "SEXCHR=${SEXCHR}" > sexchr.sh
[ "$GENDER" != "" ] && echo "GENDER=${GENDER}" > gender.sh
echo "PLATFORM=${PLATFORM}" > platform.sh

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

if [ "$splitReadArray" != "" ]; then\
	if echo $splitReadArray | grep 1 1>/dev/null;
	then
		if [ ! -e $READ1 ]
		then
			echo "FAIL Read 1 $READ1 file does not exist!"
			exit 1
		fi
		read1Size=$(ls -lah $READ1 | awk '{print $5}')
	else
		read1Size=0
	fi
	
	if echo $splitReadArray | grep 2 1>/dev/null;
	then
		if [ ! -e $READ2 ]
		then
			echo "FAIL Read 1 $READ2 file does not exist!"
			exit 1
		fi
		read2Size=$(ls -lah $READ2 | awk '{print $5}')
	else
		read2Size=0
	fi
	
	#if [ ! -e $READ1 ] || [ ! -e $READ2 ]
	#then
	#	echo "Read files may not exist!"
	#	exit 1
	#fi
	
	####################
	# Get a fancy size #
	####################
	
	
	#printf "%-22s%s" "Command" "sbatch $(dispatch "RS") -J RS_${SAMPLE}_${readSize} -a $splitReadArray $SLSBIN/readsplit.sl $SAMPLE $READ1 $READ2 $PLATFORM ${MULTI_RUN}"
	#exit 0
	# Split array contains data so run the missing split function.
	DEP_RS=$(sbatch $(dispatch "RS") -J RS_${SAMPLE}_${read1Size}_${read2Size} -a $splitReadArray $SLSBIN/readsplit.sl $SAMPLE $READ1 $READ2 $PLATFORM ${MULTI_RUN} | awk '{print $4}')
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
