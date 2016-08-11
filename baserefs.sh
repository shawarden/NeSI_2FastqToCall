#!/bin/bash

###########
# Base Reference file contains everything I don't want to retype in every script.
#
# include: "source /path/to/baserefs.sh" at the top of pretty much every script to import these values.
###########

# Nerge exit codes from piped commands so any command that fails carries to $?
set -o pipefail

# General References.
export      PROJECT=/projects/uoo00032
export    RESOURCES=${PROJECT}/Resources
export    PLATFORMS=${RESOURCES}/Capture_Platforms/GRCh37
export          BIN=${RESOURCES}/bin
export         PBIN=${BIN}/NeSI_2FastqToCall
export       SLSBIN=${PBIN}/slurm-scripts
export DESCRIPTIONS=${RESOURCES}/FastQdescriptions.txt
export       FASTQS=${PROJECT}/fastqfiles
export       BUNDLE=${RESOURCES}/broad_bundle_b37_v2.5
export        DBSNP=${BUNDLE}/dbsnp_141.GRCh37.vcf
export        MILLS=${BUNDLE}/Mills_and_1000G_gold_standard.indels.b37.vcf
export       INDELS=${BUNDLE}/1000G_phase1.indels.b37.vcf
export       COMMON=${RESOURCES}/Hapmap3_3commonvariants.vcf
export          REF=${BUNDLE}/human_g1k_v37
export         REFA=${REF}.fasta
export         REFD=${REF}.dict
export JOB_TEMP_DIR=$([ "${TMPDIR}" != "" ] && echo "${TMPDIR}" || echo "${PROJECT}/.tmp")	#system defined temp dir. /tmp/jobs/$SLURM_JOB_USER/$SLURM_JOB_ID

export      BWA=${BIN}/bwa-0.7.15/bwa
export   PICARD=${BIN}/picard-tools-2.5.0/picard.jar
export SAMTOOLS=${BIN}/samtools-1.3.1/samtools
export     GATK=${BIN}/GenomeAnalysisTK-nightly-2016-07-22-g8923599/GenomeAnalysisTK.jar

####################
# Modules versions #
####################

export MOD_ZLIB="zlib"
export MOD_JAVA="Java"

###########
# Contigs #
###########

export CONTIGS=$(cat ${REFD} | awk -F'[[:space:]:]' 'NR!=1{print $3}')		# List of contigs.
export CONTIGA=("" $(cat ${REFD} | awk -F'[[:space:]:]' 'NR!=1{print $3}'))	# Array of contigs. "" at start for 1-indexed.
export NUMCONTIGS=$(($(echo ${#CONTIGA[@]}) - 1))							# Number of contigs.

#########################
# GrCh37 Gender regions #
#########################

export XPAR1="X:1-2699520"
export TRUEX="X:2699521-154931043"
export XPAR2="X:154931044-155260560"
export TRUEY="Y:2649521-59034050"

###############################
# FastQ Split Data management #
###############################

export FASTQ_MAXREAD=15000000	# How many reads per block.
export FASTQ_MAXSCAN=10000		# How many lines to check for best index.
export FASTQ_MAXDIFF=2			# Maximum index variation before new index is created.
export FASTQ_MAXJOBS=1000		# Maximum number of alignment & sort array elements.
export FASTQ_MAXJOBZ=$(($FASTQ_MAXJOBS - 1))	# Maximum number of alignment & sort array elements starting from 0.
export FASTQ_MAXZPAD=${#FASTQ_MAXJOBS}	# Number of characters to pad to blocks.

export MAX_JOB_RATE=30			# Minimum number of seconds between job submissions.


##########################
# Dispatch function data #
##########################
declare -A SB
SB[ACCOUNT]=uoo00032
SB[MAILUSER]=sam.hawarden@otago.ac.nz
SB[MAILTYPE]=FAIL

# MWT: Max Wall-Times
# MPC: Memory per Cores
# CPT: Cores per Task

SB[RS,MWT]=0-02:30:00
SB[RS,MPC]=2048
SB[RS,CPT]=8
 
SB[PA,MWT]=0-01:30:00
SB[PA,MPC]=2048
SB[PA,CPT]=8

SB[SS,MWT]=0-01:30:00
SB[SS,MPC]=4096
SB[SS,CPT]=4

SB[CS,MWT]=0-00:15:00
SB[CS,MPC]=2048
SB[CS,CPT]=2

SB[MC,MWT]=0-02:00:00
SB[MC,MPC]=4096
SB[MC,CPT]=4

SB[MD,MWT]=0-03:00:00
SB[MD,MPC]=8192
SB[MD,CPT]=2

SB[BR,MWT]=0-01:00:00
SB[BR,MPC]=4096
SB[BR,CPT]=4

SB[PR,MWT]=0-02:00:00
SB[PR,MPC]=4092
SB[PR,CPT]=8

SB[DC,MWT]=0-00:30:00
SB[DC,MPC]=2048
SB[DC,CPT]=4

SB[GD,MWT]=0-00:10:00
SB[GD,MPC]=512
SB[GD,CPT]=1

SB[HC,MWT]=0-01:30:00
SB[HC,MPC]=4096
SB[HC,CPT]=8

SB[CR,MWT]=0-01:00:00
SB[CR,MPC]=4096
SB[CR,CPT]=1

SB[RI,MWT]=0-01:30:00
SB[RI,MPC]=4096
SB[RI,CPT]=1

SB[CV,MWT]=0-02:00:00
SB[CV,MPC]=16384
SB[CV,CPT]=1

SB[FP,MWT]=0-06:00:00
SB[FP,MPC]=4096
SB[FP,CPT]=8

SB[SV,MWT]=0-06:00:00
SB[SV,MPC]=4096
SB[SV,CPT]=8

SB[TF,MWT]=0-00:20:00
SB[TF,MPC]=512
SB[TF,CPT]=1

export SB

###########################
# Java memory calculation #
###########################

# Sometimes we call this file outside of a slurm script so fake it til it's made.
[ -z $SLURM_MEM_PER_CPU ] && SLURM_MEM_PER_CPU=4096
[ -z $SLURM_JOB_CPUS_PER_NODE ] && SLURM_JOB_CPUS_PER_NODE=4

export JAVA_MEM_GB=$(((($SLURM_MEM_PER_CPU * $SLURM_JOB_CPUS_PER_NODE)/1024)-2))
export   JAVA_ARGS="-Xmx${JAVA_MEM_GB}g -Djava.io.tmpdir=${JOB_TEMP_DIR}"
export MAX_RECORDS=$((${JAVA_MEM_GB} * 200000)) #~100bp picard records in memory.

export OPENBLAS_MAIN_FREE=1

####################
# Helper functions #
####################

###################
# Output HH:MM:SS format for a number of seconds.
###################
function printHMS {
	SECS=${1}
	printf "%02d:%02d:%02d" "$(($SECS / 3600))" "$((($SECS % 3600) / 60))" "$(($SECS % 60))"
}
export -f printHMS

#######################
# Output basic node information on job failure.
#######################
function scriptFailed {
	echo ""
	echo "$HEADER: SControl -----"
	scontrol show job ${SLURM_JOBID}
	echo ""
	echo "$HEADER: Export -----"
	export
	echo ""
	echo "$HEADER: Storage -----"
	df -ah
	echo ""
}
export -f scriptFailed

#######################
# Output runtime metrics to a log file.
#######################
function storeMetrics {
	if [ "${1}" != "" ]; then
		BACKDIR="../"
	else
		BACKDIR=""
	fi
	
	printf \
		"%19s: Task %-20s running on %2d@%1.1fGHz w/%3dGB took %s\n" \
		"$(date '+%Y-%m-%d %H:%M:%S')" \
		"${SLURM_JOB_NAME}$([ "$SLURM_ARRAY_TASK_ID" != "" ] && echo -ne ":$SLURM_ARRAY_TASK_ID")" \
		${SLURM_JOB_CPUS_PER_NODE} \
		$(echo "$(lscpu | grep "CPU MHz" | awk '{print $3}') / 1000" | bc -l) \
		$(((${SLURM_JOB_CPUS_PER_NODE} * ${SLURM_MEM_PER_CPU}) / 1024)) \
		$(printHMS $SECONDS) | \
		tee -a ${BACKDIR}metrics.txt >> ${HOME}/metrics.txt
}
export -f storeMetrics

#####################
# Return a string with a new item on the end
# Items are separated by a char or space
#####################
function appendList {
	oldList="${1}"
	newItem="${2}"
	itemJoiner=$([ "${3}" == "" ] && echo -ne " " || echo -ne "${3}")	# If blank, use space.
	
#	(>&2 echo "${oldList}${itemJoiner}${newItem}")
	
	if [ "$newItem" != "" ]; then
		if [ "$oldList" == "" ]; then
			# Initial Entry
			printf "%s" "$newItem"
		else
			# Additional Entry
			printf "%s%s%s" "$oldList" "$itemJoiner" "$newItem"
		fi
	else
		printf "%s" "$oldList"
	fi
}
export -f appendList

######################
# Return the string with the split char replaced by a space
######################
function splitByChar {
	input=${1}
	char=${2}
	
	if [ "$char" != "" ]; then
		echo -ne "$input" | sed -e "s/${char}/ /g"
	else
		echo -ne "FAIL"
		(>&2 echo -ne "FAILURE:\t splitByChar [${1}] [${2}]\n\tNo character to replace. You forget to quote your input?\n")
	fi
}
export -f splitByChar

######################
# Find matching task elements in a parent and child array.
# Set the child array task element to be dependent on the matching parent task element.
######################
function tieTaskDeps {
	childArray=${1}
	childJobID=${2}
	parentArray=${3}
	parentJobID=${4}
	
	if [ "$childArray" != "" ] && [ "$parentArray" != "" ]; then
		# Both arrays contain something.
		
		# Cycle through child array for elements
		for i in $(splitByChar "$childArray" ","); do
			# Cycle through parent array for matching elements.
			for j in $(splitByChar "$parentArray" ","); do
				if [ "$i" == "$j" ]; then
					# Match found.
#					printf " T[%s->%s] " "${childJobID}_$i" "${parentJobID}_$j"
					scontrol update JobId=${childJobID}_$i Dependency=afterok:${parentJobID}_$j
				fi
			done
		done
	fi
}
export -f tieTaskDeps

##################
# Create output olders in the temp and final destination locations.
##################
function outDirs {
	# Make sure destination location exists.
	if ! mkdir -p $(dirname ${OUTPUT}); then
		echo "$HEADER: Unable to create output folder ${PWD}/${OUTPUT}!"
		exit 1
	fi
	
	if ! mkdir -p $(dirname ${JOB_TEMP_DIR}/${OUTPUT}); then
		echo "$HEADER: Unable to create temp output folder ${JOB_TEMP_DIR}/${OUTPUT}!"
		exit 1
	fi
}
export -f outDirs

##################
# Check if final output exists.
##################
function outFile {
	if [ -e ${OUTPUT} ]; then
		echo "$HEADER: Output file \"${OUTPUT}\" already exists!"
		exit 1
	fi
}
export -f outFile

#################
# Make sure input file exists
#################
function inFile {
	if [ ! -e ${INPUT} ]; then
		echo "$HEADER: Input file \"${INPUT}\" doesn't exists!"
		ls -la $(dirname $INPUT)
		exit 1
	fi
}
export inFile

###################
# Move temp output to final output locations
###################
function finalOut {
	if ! mv ${JOB_TEMP_DIR}/${OUTPUT%.*}* $(dirname ${OUTPUT}); then
		echo "$HEADER: Failed to move ${JOB_TEMP_DIR}/${OUTPUT} to ${PWD}/${OUTPUT}!"
		exit 1
	fi
}
export -f finalOut

####################
# Check if command filed
####################
function cmdFailed {
	echo "$HEADER: [$SLURM_JOB_NAME:$SLURM_JOBID:$SLURM_ARRAY_TASK_ID] failed with $?!"
}
export -f cmdFailed

###################
# Outputs a dependency if one exists.
###################
function depCheck {
	if [ "$1" != "" ]
	then
		echo -ne "--dependency afterok:${1}"
	fi
}
export -f depCheck

###################
# Gets job category data
###################
function dispatch {
	jobWait
	JOBCAT=${1}
	echo -ne "--account ${SB[ACCOUNT]} --mail-user=${SB[MAILUSER]} --mail-type=${SB[MAILTYPE]} --time ${SB[${JOBCAT},MWT]} --mem-per-cpu=${SB[${JOBCAT},MPC]} --cpus-per-task ${SB[${JOBCAT},CPT]}"
}
export -f dispatch

###################
# Echo job stats to log file
###################
function jobStats {
	echo -e "$HEADER:\tMax Walltime:     ${SB[$HEADER,MWT]}"
	echo -e "\tMemory-Per-Core: ${SB[$HEADER,MPC]}"
	echo -e "\tCores-Per-Task:  ${SB[$HEADER,CPT]}"
}
export -f jobStats

##################
# Delays job submission based on project's previous submitted job.
##################
function jobWait {
	timeNow=$(date +%s)
	if [ -e ${PROJECT}/submitted.log ]; then
		lastJob=$(cat ${PROJECT}/submitted.log)
		if [ "$lastJob" != "" ]; then
			>&2 printf " "
			while [ $(($timeNow - $lastJob)) -lt $MAX_JOB_RATE ]; do
				case "$TICKER" in
					"|") TICKER="/" ;;
					"/") TICKER="-" ;;
					"-") TICKER="\\" ;;
					*) TICKER="|" ;;
				esac
				
				>&2 printf "\b%s" "$TICKER"
				
				sleep 0.25s
				timeNow=$(date +%s)
			done
			>&2 printf "\b"
		fi
	fi
	echo $timeNow > ${PROJECT}/submitted.log
}
