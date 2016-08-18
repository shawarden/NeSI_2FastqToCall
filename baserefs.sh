#!/bin/bash

###########
# Base Reference file contains everything I don't want to retype in every script.
#
# include: "source /path/to/baserefs.sh" at the top of pretty much every script to import these values.
###########

# Nerge exit codes from piped commands so any command that fails carries to $?
set -o pipefail

##############
# Exit codes #
##############

export EXIT_IO=10
export EXIT_PR=15
export EXIT_MV=20
export EXIT_TF=21

######################
# General References #
######################

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
export JOB_TEMP_DIR=$([ "${TMPDIR}" != "" ] && echo "${TMPDIR}" || echo "${PROJECT}/.tmp")

##################################
# Version controlled executables #
##################################

export      BWA=${BIN}/bwa-0.7.15/bwa
export   PICARD=${BIN}/picard-tools-2.5.0/picard.jar
export SAMTOOLS=${BIN}/samtools-1.3.1/samtools
export     GATK=${BIN}/GenomeAnalysisTK-nightly-2016-07-22-g8923599/GenomeAnalysisTK.jar

######################
# Contig shinanigans #
######################

# Setup basic contig definition.
export        REFD=${REF}.dict
export CONTIGARRAY=("" $(cat ${REFD} | awk -F'[[:space:]:]' 'NR!=1{print $3}'))
export  NUMCONTIGS=$(($(echo ${#CONTIGARRAY[@]}) - 1))

# Gender vs autosomal contigs.
export  SEXCONTIGS="23,24"	# Gender contigs that are special cases.
export AUTOCONTIGS="1-22"	# Autosomal contigs that determine base coverage for gender calculation.

# Gender contigs have parts that are not created equal.
export       XPAR1="X:1-2699520"
export       TRUEX="X:2699521-154931043"
export       XPAR2="X:154931044-155260560"
export       TRUEY="Y:2649521-59034050"


####################
# Modules versions #
####################

export MOD_ZLIB="zlib"
export MOD_JAVA="Java"

###############################
# FastQ Split Data management #
###############################

export FASTQ_MAXREAD=15000000	# How many reads per block.
export FASTQ_MAXSCAN=10000		# How many lines to check for best index.
export FASTQ_MAXDIFF=2			# Maximum index variation before new index is created.
export FASTQ_MAXJOBS=100		# Maximum number of alignment & sort array elements.
export FASTQ_MAXJOBZ=99	#$(($FASTQ_MAXJOBS - 1))	# Maximum number of alignment & sort array elements starting from 0.
export FASTQ_MAXZPAD=4	#${#FASTQ_MAXJOBS}	# Number of characters to pad to blocks.

# Minimum number of seconds between job submissions
# 200 submissions/10minutes = 20 submissions/minute = 3 seconds/submission.
# 600 submissions/60minutes = 10 submissions/minute = 6 seconds/submission.
# 5000 total jobs (including all array elements.
export MAX_JOB_RATE=6


##########################
# Dispatch function data #
##########################
declare -A SB
SB[ACCOUNT]=uoo00032
SB[MAILUSER]=sam.hawarden@otago.ac.nz
SB[MAILTYPE]="ARRAY_TASKS,TIME_LIMIT_90,FAIL"

# MWT: Calibrated Max Wall-Time
# MPC: Memory per Cores
# CPT: Cores per Task

# Calculating base and per contig walltimes.
# combined read filesize:
# WALLTIME=BASE_MULTIPLIER * CONTIG_MULTIPLIER * CALIBRATED_MAX_WALLTIME * (INPUT_SIZE / CALIBRATED_INPUT_SIZE)
# WALLTIME=$(printf "%f\n" $(echo "${SB[BWTM]} * ${SB[WTM,${CONTIGNUMBER}]} * ${SB[${HEADER},MWT]} * ($INPUT_FILE_SIZE / $CALIBRATED_FILE_SIZE)" | bc -l))
#

SB[BWTM]=1.25	# Base wall-time multiplier.

# ReadSplit.
SB[RS,MWT]=210
SB[RS,MPC]=2048
SB[RS,CPT]=8

# PrimaryAlignment.
SB[PA,MWT]=90
SB[PA,MPC]=2048
SB[PA,CPT]=8

# SortSAM
SB[SS,MWT]=90
SB[SS,MPC]=4096
SB[SS,CPT]=4

# ContigSplit
SB[CS,MWT]=30
SB[CS,MPC]=2048
SB[CS,CPT]=1

# MergeContigs
SB[MC,MWT]=120
SB[MC,MPC]=4096
SB[MC,CPT]=4

# MarkDuplicates
SB[MD,MWT]=180
SB[MD,MPC]=8192
SB[MD,CPT]=2

# BaseRecalibration
SB[BR,MWT]=60
SB[BR,MPC]=4096
SB[BR,CPT]=4

# PrintReads
SB[PR,MWT]=120
SB[PR,MPC]=4092
SB[PR,CPT]=8

# DepthofCoverage
SB[DC,MWT]=60
SB[DC,MPC]=2048
SB[DC,CPT]=8

# GenderDetermination
SB[GD,MWT]=10
SB[GD,MPC]=512
SB[GD,CPT]=1

# HaplotypeCaller
SB[HC,MWT]=90
SB[HC,MPC]=4096
SB[HC,CPT]=8

# CatReads
SB[CR,MWT]=60
SB[CR,MPC]=4096
SB[CR,CPT]=1

# ReadIndex
SB[RI,MWT]=93
SB[RI,MPC]=4096
SB[RI,CPT]=1

# CatVariants
SB[CV,MWT]=120
SB[CV,MPC]=16384
SB[CV,CPT]=1

# FingerPrint
SB[FP,MWT]=360
SB[FP,MPC]=4096
SB[FP,CPT]=8

# SelectVariants
SB[SV,MWT]=360
SB[SV,MPC]=4096
SB[SV,CPT]=8

# TransferFile
SB[TF,MWT]=20
SB[TF,MPC]=512
SB[TF,CPT]=1

#Wall time multipliers to dynamically decrease array element walltime.
SB[WTM,1]=1.0
SB[WTM,2]=${SB[WTM,1]}
SB[WTM,3]=1.0
SB[WTM,4]=1.0
SB[WTM,5]=1.0
SB[WTM,6]=1.0
SB[WTM,7]=1.0
SB[WTM,8]=1.0
SB[WTM,9]=1.0
SB[WTM,10]=1.0
SB[WTM,11]=1.0
SB[WTM,12]=1.0
SB[WTM,13]=1.0
SB[WTM,14]=1.0
SB[WTM,15]=1.0
SB[WTM,16]=1.0
SB[WTM,17]=1.0
SB[WTM,18]=1.0
SB[WTM,19]=1.0
SB[WTM,20]=1.0
SB[WTM,21]=1.0
SB[WTM,22]=1.0
SB[WTM,23]=1.0
SB[WTM,24]=1.0
SB[WTM,25]=1.0
SB[WTM,26]=${SB[WTM,26]}
SB[WTM,27]=${SB[WTM,26]}
SB[WTM,28]=${SB[WTM,26]}
SB[WTM,29]=${SB[WTM,26]}
SB[WTM,30]=${SB[WTM,26]}
SB[WTM,31]=${SB[WTM,26]}
SB[WTM,32]=${SB[WTM,26]}
SB[WTM,33]=${SB[WTM,26]}
SB[WTM,34]=${SB[WTM,26]}
SB[WTM,35]=${SB[WTM,26]}
SB[WTM,36]=${SB[WTM,26]}
SB[WTM,37]=${SB[WTM,26]}
SB[WTM,38]=${SB[WTM,26]}
SB[WTM,39]=${SB[WTM,26]}
SB[WTM,40]=${SB[WTM,26]}
SB[WTM,41]=${SB[WTM,26]}
SB[WTM,42]=${SB[WTM,26]}
SB[WTM,43]=${SB[WTM,26]}
SB[WTM,44]=${SB[WTM,26]}
SB[WTM,45]=${SB[WTM,26]}
SB[WTM,46]=${SB[WTM,26]}
SB[WTM,47]=${SB[WTM,26]}
SB[WTM,48]=${SB[WTM,26]}
SB[WTM,49]=${SB[WTM,26]}
SB[WTM,50]=${SB[WTM,26]}
SB[WTM,51]=${SB[WTM,26]}
SB[WTM,52]=${SB[WTM,26]}
SB[WTM,53]=${SB[WTM,26]}
SB[WTM,54]=${SB[WTM,26]}
SB[WTM,55]=${SB[WTM,26]}
SB[WTM,56]=${SB[WTM,26]}
SB[WTM,57]=${SB[WTM,26]}
SB[WTM,58]=${SB[WTM,26]}
SB[WTM,59]=${SB[WTM,26]}
SB[WTM,60]=${SB[WTM,26]}
SB[WTM,61]=${SB[WTM,26]}
SB[WTM,62]=${SB[WTM,26]}
SB[WTM,63]=${SB[WTM,26]}
SB[WTM,64]=${SB[WTM,26]}
SB[WTM,65]=${SB[WTM,26]}
SB[WTM,66]=${SB[WTM,26]}
SB[WTM,67]=${SB[WTM,26]}
SB[WTM,68]=${SB[WTM,26]}
SB[WTM,69]=${SB[WTM,26]}
SB[WTM,70]=${SB[WTM,26]}
SB[WTM,71]=${SB[WTM,26]}
SB[WTM,72]=${SB[WTM,26]}
SB[WTM,73]=${SB[WTM,26]}
SB[WTM,74]=${SB[WTM,26]}
SB[WTM,75]=${SB[WTM,26]}
SB[WTM,76]=${SB[WTM,26]}
SB[WTM,77]=${SB[WTM,26]}
SB[WTM,78]=${SB[WTM,26]}
SB[WTM,79]=${SB[WTM,26]}
SB[WTM,80]=${SB[WTM,26]}
SB[WTM,81]=${SB[WTM,26]}
SB[WTM,82]=${SB[WTM,26]}
SB[WTM,83]=${SB[WTM,26]}
SB[WTM,84]=${SB[WTM,26]}

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

##
# Outputs minutes from various combinations of time strings:
#  D-HH:MM:SS
#    HH:MM:SS
#       MM:SS
#       MM
##
function printMinutes {
	INPUT="${1}"
	DAYS=$(echo ${INPUT} | cut -d'-' -f1)
	[ "$DAYS" == "" ] && DAYS=0
	
	INPUT=$(echo ${INPUT} | cut -d'-' -f2)
	local numTimeBlocks=$(($(echo $INPUT | grep -o ':' | wc -l) + 1))
	case $numTimeBlocks in
		1)	# Minutes
			MINS=$INPUT
			;;
		2)	# Minutes:Seconds
			MINS=$(echo $INPUT | cut -d':' -f1)
			SECS=$(echo $INPUT | cut -d':' -f2)
			;;
		3)	# Hours:Minutes:Seconds
			HOURS=$(echo $INPUT | cut -d':' -f1)
			 MINS=$(echo $INPUT | cut -d':' -f2)
			 SECS=$(echo $INPUT | cut -d':' -f3)
			;;
		*)
			echo "WHAT?"
			exit 1
	esac
	
	printf "%d" $((($DAYS * 1440) + ($HOURS * 60) + $MINS + ($SECS/60)))
}
export -f printMinutes

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
function failMetrics {
	SIGTERM="SIGTERM"
	
#	case $HEADER in
#		RS|PA|SS|CS)
#			BACKDIR="../"
#			;;
#		*)
#			BACKDIR=""
#	esac
#	
#	printf \
#		"%s\t%s\t%dc\t%1.1fGHz\t%dGB\t%s\tSIGTERM\n" \
#		"$(date '+%Y-%m-%d %H:%M:%S')" \
#		"${SLURM_JOB_NAME}$([ "$SLURM_ARRAY_TASK_ID" != "" ] && echo -ne ":$SLURM_ARRAY_TASK_ID")" \
#		${SLURM_JOB_CPUS_PER_NODE} \
#		$(echo "$(lscpu | grep "CPU MHz" | awk '{print $3}') / 1000" | bc -l) \
#		$(((${SLURM_JOB_CPUS_PER_NODE} * ${SLURM_MEM_PER_CPU}) / 1024)) \
#		$(printHMS $SECONDS) | \
#		tee -a ${BACKDIR}metrics.txt | \
#		tee -a ${HOME}/metrics.txt >> ${HOME}/$(date '+%Y_%m_%d').metrics.txt
}
export -f failMetrics

trap "failMetrics" SIGTERM

#######################
# Output runtime metrics to a log file.
#######################
function storeMetrics {
	if [ "$SLURM_JOB_NAME" != "" ] && [ "$HEADER" != "" ] ; then
		case $HEADER in
			RS|PA|SS|CS)
				BACKDIR="../"
				;;
			*)
				BACKDIR=""
		esac
		
		printf \
			"%s\t%s\t%dc\t%1.1fGHz\t%dGB\t%s%s\n" \
			"$(date '+%Y-%m-%d %H:%M:%S')" \
			"${SLURM_JOB_NAME}$([ "$SLURM_ARRAY_TASK_ID" != "" ] && echo -ne ":$SLURM_ARRAY_TASK_ID")" \
			${SLURM_JOB_CPUS_PER_NODE} \
			$(echo "$(lscpu | grep "CPU MHz" | awk '{print $3}') / 1000" | bc -l) \
			$(((${SLURM_JOB_CPUS_PER_NODE} * ${SLURM_MEM_PER_CPU}) / 1024)) \
			$(printHMS $SECONDS) \
			$([ "$SIGTERM" != "" ] && echo -ne "	SIGTERM" || echo -ne "") | \
			tee -a ${BACKDIR}metrics.txt | \
			tee -a ${HOME}/metrics.txt >> ${HOME}/$(date '+%Y_%m_%d').metrics.txt
	fi
}
export -f storeMetrics

trap "storeMetrics" EXIT

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
		# Cycle through child array elements
		for i in $(splitByChar "$childArray" ","); do
			elementMatched=0
			# Cycle through parent array elements.
			for j in $(splitByChar "$parentArray" ","); do
				if [ "$i" == "$j" ]; then
					# Matching element found. Tie child element to parent element.
#					printf " T[%s->%s] " "${childJobID}_$i" "${parentJobID}_$j"
					scontrol update JobId=${childJobID}_$i Dependency=afterok:${parentJobID}_$j
					elementMatched=1
				fi
			done
			if [ $elementMatched -eq 0 ]; then
				# No matching element found in parent array.
				# Release child element from entire parent array.
				scontrol update JobId=${childJobID}_$i Dependency=
			fi
			tickOver
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
	SIGTERM="FAIL"
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
	printf "%s--account %s --mail-user %s --mail-type %s --time %.0f --mem-per-cpu %d --cpus-per-task %d" \
	"" \
	${SB[ACCOUNT]} \
	${SB[MAILUSER]} \
	${SB[MAILTYPE]} \
	${SB[${JOBCAT},MWT]} \
	${SB[${JOBCAT},MPC]} \
	${SB[${JOBCAT},CPT]}
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
			while [ $(($timeNow - $lastJob)) -lt $MAX_JOB_RATE ]; do
				tickOver
				
				sleep 0.25s
				timeNow=$(date +%s)
			done
		fi
	fi
	echo $timeNow > ${PROJECT}/submitted.log
}

###################
# Visual ticker
###################
function tickOver {
	case "$TICKER" in
		"|") TICKER="/" ;;
		"/") TICKER="-" ;;
		"-") TICKER="\\" ;;
		*) TICKER="|" ;;
	esac
	
	>&2 printf "%s\b" "$TICKER"
}
export -f tickOver
