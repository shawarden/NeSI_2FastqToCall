#!/bin/bash

###########
# Base Reference file contains everything I don't want to retype in every script.
#
# include: "source /path/to/baserefs.sh" at the top of pretty much every script to import these values.
###########

# Merge exit codes from piped commands so any command that fails carries to $?
set -o pipefail

##############
# Exit codes #
##############

export EXIT_IO=12
export EXIT_PR=15
export EXIT_MV=20
export EXIT_TF=21

######################
# General References #
######################

export SLURM_VERSION=$(scontrol -V | awk '{print $2}')


export ENDPOINT_NESI=d62d1ead-6d04-11e5-ba46-22000b92c6ec
export ENDPOINT_UOO=7770a464-0540-11e6-a72d-22000bf2d559

export      PROJECT=/projects/uoo00032
export    RESOURCES=${PROJECT}/Resources
export    PLATFORMS=${RESOURCES}/Capture_Platforms/GRCh37
export          BIN=${RESOURCES}/bin
export         PBIN=${BIN}/NeSI_2FastqToCall
export       SLSBIN=${PBIN}/slurm-scripts
export DESCRIPTIONS=${RESOURCES}/FastQdescriptions.txt
export       FASTQS=${PROJECT}/fastqfiles

#export       COMMON=${RESOURCES}/Hapmap3_3commonvariants.vcf
export       BUNDLE=${RESOURCES}/broad_bundle_b37_v2.5
export        DBSNP=${BUNDLE}/dbsnp_141.GRCh37.vcf
export        MILLS=${BUNDLE}/Mills_and_1000G_gold_standard.indels.b37.vcf
export       INDELS=${BUNDLE}/1000G_phase1.indels.b37.vcf
export          REF=${BUNDLE}/human_g1k_v37_decoy

export       COMMON=${RESOURCES}/Hapmap3_3commonvariants.vcf
#export       BUNDLE=${RESOURCES}/v0
#export        DBSNP=${BUNDLE}/Homo_sapiens_assembly38.dbsnp138.vcf
#export        MILLS=${BUNDLE}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#export       INDELS=${BUNDLE}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
#export          REF=${BUNDLE}/Homo_sapiens_assembly38

export         REFA=${REF}.fasta
export JOB_TEMP_DIR=$([ "${TMPDIR}" != "" ] && echo "${TMPDIR}" || echo "${PROJECT}/.tmp")

##################################
# Version controlled executables #
##################################

export       BWA=${BIN}/bwa-0.7.15/bwa
export    PICARD=${BIN}/picard-tools-2.5.0/picard.jar
export  SAMTOOLS=${BIN}/samtools-1.3.1/samtools
#export      GATK=${BIN}/GenomeAnalysisTK-nightly-2016-07-22-g8923599/GenomeAnalysisTK.jar
export      GATK=${BIN}/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
export   ZIP_CMD=${BIN}/pigz-2.3.3/pigz_sb
export   CAT_CMD="${ZIP_CMD} -cd"
export SPLIT_CMD=${BIN}/coreutils-8.25/src/split

######################
# Contig shinanigans #
######################

# Setup basic contig definition.
export REFD=${REF}.dict

if [ ! -e $REFD ]; then 
	echo "Unable to locate reference dictionary"
	exit 1
fi

# List of contigs in specified reference sequence.
# GRCh37 has 84
# HG38 has 3367
export CONTIGARRAY=("" $(cat $REFD | awk 'NR!=1{print $2}' | sed -e 's/SN://g'))

# List of contigs in specified reference sequence, excluding those marked as decoy, non-human sequences.
# GRCh37 has 84
# HG38 has 982 (2384 decoys removed)
export CONTIGARRAY_NODECOY=("" $(cat $REFD | grep -v "decoy" | awk 'NR!=1{print $2}' | sed -e 's/SN://g'))

# List of primary contigs to wrap all alternate and positionless contigs within a chromosome.
# In HG38 chr1 will be combination of chr1, chr1_KI270706v1_random, chr1_KI270762v1_alt, etc...
# In HG38 HLA will contain all HLA contigs.
# In HG38 EBV will contain all Epstein-Barr Virus contigs
# GRCh37 has 84
# HG38 has 28
export CONTIGBLOCKS=("" $(for i in $(seq 0 $((${#CONTIGARRAY[@]}-1))); do echo "${CONTIGARRAY[$i]}"; done | awk -F'[* ]' '
{
	words[$1]++
}
END{
	for (w in words) {
		printf "%s\n", w
	}
}
	' | sort -V -f
))

export         NUMCONTIGS=$((${#CONTIGARRAY[@]} - 1))
export NUMCONTIGS_NODECOY=$((${#CONTIGARRAY_NODECOY[@]} - 1))
export   NUMCONTIG_BLOCKS=$((${#CONTIGBLOCKS[@]} - 1))

# Gender contigs have parts that are not created equal.
export XPAR1="X:1-2699520"
export TRUEX="X:2699521-154931043"
export XPAR2="X:154931044-155260560"
export TRUEY="Y:2649521-59034050"

#export YPAR1="chrY:10000-2781479"		# Hard masked to Ns in HG38
#export YPAR2="chrY:56887902-57217415"
#export TRUEY="chrY:2781480-56887901"	# Hard masked to Ns in HG38

#export XPAR1="chrX:10000-2781479"
#export XPAR2="chrX:155701382-156030895"
#export TRUEX="chrX:2781480-155701381"

####################
# Modules versions #
####################

export MOD_ZLIB="zlib"
export MOD_JAVA="Java"

###############################
# FastQ Split Data management #
###############################

export FASTQ_MAXREAD=10000000	# How many reads per block.
export FASTQ_MAXSCAN=10000		# How many lines to check for best index.
export FASTQ_MAXDIFF=2			# Maximum index variation before new index is created.
export FASTQ_MAXJOBS=100		# Maximum number of alignment & sort array elements.
export FASTQ_MAXJOBZ=99			#$(($FASTQ_MAXJOBS - 1))	# Maximum number of alignment & sort array elements starting from 0.
export FASTQ_MAXZPAD=4			#${#FASTQ_MAXJOBS}	# Number of characters to pad to blocks.

# Minimum number of seconds between job submissions
# 200 submissions/10minutes = 20 submissions/minute = 3 seconds/submission.
# 600 submissions/60minutes = 10 submissions/minute = 6 seconds/submission.
# 5000 total jobs (including all array elements.
# Set to zero as merged align, sort and split job into one so no submission rate issue.
export MAX_JOB_RATE=6


##########################
# Dispatch function data #
##########################
declare -A SB
SB[ACCOUNT]=uoo00032
SB[MAILUSER]=david.markie@otago.ac.nz
SB[MAILTYPE]=FAIL

# MWT: Calibrated Max Wall-Time
# MPC: Memory per Cores
# CPT: Cores per Task

# Calculating base and per contig walltimes.
# combined read filesize:
# WALLTIME=BASE_MULTIPLIER * CONTIG_MULTIPLIER * CALIBRATED_MAX_WALLTIME * (INPUT_SIZE / CALIBRATED_INPUT_SIZE)
# WALLTIME=$(printf "%f\n" $(echo "${SB[BWTM]} * ${SB[WTM,${CONTIGNUMBER}]} * ${SB[${HEADER},MWT]} * ($INPUT_FILE_SIZE / $CALIBRATED_FILE_SIZE)" | bc -l))
#

SB[BWTM]=1.25	# Base wall-time multiplier.
SB[MWT]=359		# Maximum wall-time to be in High partition.
#SB[MWT]=1440	# 1 day.

# ReadSplit.
SB[RS]="RSplit"
SB[RS,MWT]=${SB[MWT]}
SB[RS,MPC]=512
SB[RS,CPT]=6

# Block alignment (Merged PA, SS & CS)
SB[PA]="Align"
SB[SS]="Sort"
SB[CS]="CSplit"
SB[BA]="AlignSortSplit"
SB[BA,MWT]=${SB[MWT]}
SB[BA,MPC]=2048
SB[BA,CPT]=8

# Merge & Mark
SB[MC]="CMerge"
SB[MD]="MarkDup"
SB[MM]="MergeMark"
SB[MM,MWT]=719
SB[MM,MPC]=8192
SB[MM,CPT]=2

# Recal & Print
SB[BR]="BaseRecal"
SB[PR]="PrintReads"
SB[RC]="ReCal"
SB[RC,MWT]=${SB[MWT]}
SB[RC,MPC]=4092
SB[RC,CPT]=8

# DepthofCoverage
SB[DC]="Depth"
SB[DC,MWT]=${SB[MWT]}
SB[DC,MPC]=2048
SB[DC,CPT]=8

# GenderDetermination
SB[GD]="Gender"
SB[GD,MWT]=${SB[MWT]}
SB[GD,MPC]=512
SB[GD,CPT]=1

# HaplotypeCaller
SB[HC]="Haplotype"
SB[HC,MWT]=${SB[MWT]}
SB[HC,MPC]=4096
SB[HC,CPT]=8

# CatReads
SB[CR]="CatReads"
SB[CR,MWT]=${SB[MWT]}
SB[CR,MPC]=4096
SB[CR,CPT]=1

# ReadIndex
SB[RI]="IndexReads"
SB[RI,MWT]=${SB[MWT]}
SB[RI,MPC]=4096
SB[RI,CPT]=2

# CatVariants
SB[CV]="CatVar"
SB[CV,MWT]=${SB[MWT]}
SB[CV,MPC]=8192
SB[CV,CPT]=2

# FingerPrint
SB[FP]="FingerPrint"
SB[FP,MWT]=${SB[MWT]}
SB[FP,MPC]=4096
SB[FP,CPT]=8

# SelectVariants
SB[SV]="SelectVar"
SB[SV,MWT]=${SB[MWT]}
SB[SV,MPC]=4096
SB[SV,CPT]=8

# TransferFile
SB[TF]="Transfer"
SB[TF,MWT]=${SB[MWT]}
SB[TF,MPC]=1024
SB[TF,CPT]=6

export SB

###########################
# Java memory calculation #
###########################

# Sometimes we call this file outside of a slurm script so fake it til it's made.
[ -z $SLURM_MEM_PER_CPU ] && SLURM_MEM_PER_CPU=4096
[ -z $SLURM_JOB_CPUS_PER_NODE ] && SLURM_JOB_CPUS_PER_NODE=4
[ -z $RAMDISK ] && RAMDISK=0

export JAVA_MEM_GB=$(((($SLURM_MEM_PER_CPU * $SLURM_JOB_CPUS_PER_NODE) / 1024) - 3 - $RAMDISK))
export   JAVA_ARGS="-Xmx${JAVA_MEM_GB}g -Djava.io.tmpdir=${JOB_TEMP_DIR}"
export MAX_RECORDS=$((${JAVA_MEM_GB} * 200000)) #~100bp picard records in memory.

export OPENBLAS_MAIN_FREE=1

####################
# Picard Arguments #
####################

# Shared Picard arguments.
PIC_ARGS="CREATE_INDEX=true \
COMPRESSION_LEVEL=9 \
MAX_RECORDS_IN_RAM=${MAX_RECORDS} \
TMP_DIR=${JOB_TEMP_DIR}"

# Sort specific arguments
SORT_ARGS="SORT_ORDER=coordinate"

# Merge specific arguments
MERGE_ARGS="USE_THREADING=true"

# MarkDuplicate specific arguments.
MARK_ARGS="METRICS_FILE=metrics.txt"

##################
# GATK arguments #
##################

GATK_ARGS="-R ${REFA} \
-nct ${SLURM_JOB_CPUS_PER_NODE}"

GATK_BSQR="-T BaseRecalibrator \
-knownSites ${DBSNP} \
-knownSites ${MILLS}"

GATK_READ="-T PrintReads \
--bam_compression 9"

GATK_HTC="-T HaplotypeCaller \
--emitRefConfidence GVCF \
--dbsnp ${DBSNP}"

#-G StandardAnnotation
#-G AS_StandardAnnotation
#-G StandardHCAnnotation"

#--pedigree /path/to/Ped.fancydatefunctionhere.txt
#--pedigreeValidationType SILENT"

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
#  D-H:M:S
#    H:M:S
#      M:S
#      M
##
function printMinutes {
	printf "%.0f" $(echo "$(printSeconds $*) / 60" | bc -l)
}
export -f printMinutes

##
# Outputs seconds from various combinations of time strings:
#  D-H:M:S
#    H:M:S
#      M:S
#      M
##
function printSeconds {
	declare SECONDSINPUT=${*:-$(</dev/stdin)}
	echo $SECONDSINPUT | awk '{
		numHMSBlocks=split($0,hmsBlocks,":")
		
		if (numHMSBlocks == 1) {
			# Minutes only
			printf "%.0f", hmsBlocks[1] * 60
		} else if (numHMSBlocks == 2) {
			# Minutes:Seconds
			printf "%.0f", (hmsBlocks[1] * 60) + hmsBlocks[2]
		} else if (numHMSBlocks == 3) {
			# (days?)-Hours:Minutes:Seconds
			numDHBlocks=split(hmsBlocks[1],dhBlocks,"-")
			if (numDHBlocks == 1) {
				# Hours only.
				printf "%.0f", (dhBlocks[1] * 60 * 60) + (hmsBlocks[2] * 60) + hmsBlocks[3]
			} else {
				# Days-Hours.
				printf "%.0f", (dhBlocks[1] * 24 * 60 * 60) + (dhBlocks[2] * 60 * 60) + (hmsBlocks[2] * 60) + hmsBlocks[3]
			}
		}
	}'
}
export -f printSeconds

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
	export SIGTERM="SIGTERM"
	storeMetrics
}
export -f failMetrics

trap "failMetrics" SIGTERM

#######################
# Output runtime metrics to a log file.
#######################
function storeMetrics {
	sleep 10s
	if [ "$SLURM_JOB_NAME" != "" ] && [ "$HEADER" != "" ] && [ "$HEADER" != "CU" ]; then
		case $HEADER in
			RS|PA|SS|CS|BA)
				BACKDIR="../"
				;;
			*)
				BACKDIR=""
		esac
		
		if [ ! -e ${BACKDIR}metrics.txt ]; then
			echo "$logLine" > ${BACKDIR}metrics.txt
		fi
		
		logLine="DateTime,ID,RunTime,Node,CoreAlloc,CoreUtil,MemAlloc,MemUtil,Result,JobName"
		if [ ! -e $HOME/metrics.txt ]; then
			echo "$logLine" > $HOME/metrics.txt
		fi
		
		if [ ! -e $HOME/$(date '+%Y_%m_%d').metrics.txt ]; then
			echo "$logLine" > $HOME/$(date '+%Y_%m_%d').metrics.txt
		fi
		
		jobID=$([ "$SLURM_ARRAY_JOB_ID" != "" ] && echo -ne "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" || echo -ne "${SLURM_JOBID}")
		jobString="${jobID}$([ "$JOBSTEP" != "" ] && echo -ne ".${JOBSTEP}")"
		(echo "JobID[.JobStep]:$jobString" 1>&2)
		
		sacct --format jobid%20,jobname%20,elapsed,AveCPU,MinCPU,TotalCPU,UserCPU,CPUTime,CPUTimeRaw -j ${jobID}
		
		jobStats=$(sacct --format JobID%20,JobName%10,UserCPU,CPUTimeRaw,MaxRSS -j $jobString)
		(echo "$jobStats" 1>&2)
		
		JobName=${SLURM_JOB_NAME#*_}
		
		jobLine=$(echo "$jobStats" | grep "\s${jobString}\s")
		(echo "jobline: $jobLine" 1>&2)
		
		UserCPU=$(printSeconds $(echo $jobLine | awk '{print $3}'))	# Convert usercpu time to seconds.
		IdealCPU=$(echo $jobLine | awk '{print $4}')	# CPUTimeRaw is already in seconds.
		CPUUsage=$(echo "$UserCPU $IdealCPU $SLURM_JOB_CPUS_PER_NODE" | awk '{printf "%.2f", (($1 / $2) * $3)}' )
		
		(echo "CPU: ($UserCPU / $IdealCPU) * $SLURM_JOB_CPUS_PER_NODE = $CPUUsage" 1>&2)
		
		MaxRSS=$(echo $jobLine | awk '{print $5}')
		MaxRSSMB=$(echo "${MaxRSS%?}" | awk '{printf "%f", $1 / 1024}')
		MaxMem=$(echo "${SLURM_JOB_CPUS_PER_NODE} ${SLURM_MEM_PER_CPU}" | awk '{printf "%.0f", $1 * $2}' )
		MemUsage=$(echo "$MaxRSSMB $MaxMem" | awk '{printf "%.2f", ($1 / $2)}' )
		(echo "MEM: ($MaxRSSMB / $MaxMem) = $MemUsage" 1>&2)
		
		# Desired log output: DATE.Time, JobID, RunTime, CPUs, CPU Usage, MaxMem, Mem Usage, Completion state, Job Name.
		printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
			"$(date '+%Y-%m-%d %H:%M:%S')" \
			$jobString \
			$(printHMS $SECONDS) \
			${SLURM_NODELIST//compute-/} \
			$SLURM_JOB_CPUS_PER_NODE \
			$([ "$CPUUsage" != "" ] && echo -ne "$CPUUsage" || echo -ne "0.0") \
			$([ "$MaxMem" != "" ] && echo -ne "$MaxMem" || echo -ne "0.0") \
			$([ "$MaxRSSMB" != "" ] && echo -ne "$MaxRSSMB" || echo -ne "0.0") \
			$([ "$SIGTERM" != "" ] && echo -ne "$SIGTERM" || echo -ne "PASS") \
			${SB[$HEADER]} \
			"${JobName}$([ "$SLURM_ARRAY_JOB_ID" != "" ] && echo -ne ":$SLURM_ARRAY_TASK_ID")" | \
			tee -a ${BACKDIR}metrics.txt | \
			tee -a $HOME/metrics.txt >> \
			$HOME/$(date '+%Y_%m_%d').metrics.txt
	fi
}
export -f storeMetrics

#trap "storeMetrics" EXIT

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
	if [ -e ${OUTPUT}.done ]; then
		# Done file exists. why are we here?
		echo "$HEADER: Output file \"${OUTPUT}\" already completed!"
		exit 0
	elif [ -e ${OUTPUT} ]; then
		# Output already exists for this process. Overwrite!
		echo "$HEADER: Output file \"${OUTPUT}\" already exists. Overwriting!"
	fi
}
export -f outFile

#################
# Make sure input file exists
#################
function inFile {
	if [ ! -e "${INPUT}" ]; then
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

###################
# Move temp output to final output locations
###################
function scratchOut {
	if ! mv ${SCRATCH_DIR}/${OUTPUT%.*}* $(dirname ${OUTPUT}); then
		echo "$HEADER: Failed to move ${SCRATCH_DIR}/${OUTPUT} to ${PWD}/${OUTPUT}!"
		exit 1
	fi
}
export -f scratchOut

####################
# Check if command filed
####################
function cmdFailed {
	exitCode=$1
	echo "$HEADER: [$SLURM_JOB_NAME:$SLURM_JOBID:$SLURM_ARRAY_TASK_ID] failed with $exitCode!"
	SIGTERM="FAIL($exitCode)"
}
export -f cmdFailed

###################
# Outputs a dependency if one exists.
###################
function depCheck {
	(echo "depCheck: $0 \"${@}\"" >> ~/depCheck.log)
	#[ "$1" != "" ] && echo -ne "--dependency $([ "$2" != "" ] && printf "aftercorr" || printf "afterok"):${1}"
	[ "${@}a" != "a" ] && echo -ne "--dependency afterok:${@}"
}
export -f depCheck

function depCheckArrack {
	(echo "depCheck: $0 \"${@}\"" >> ~/depCheck.log)
	#[ "$1" != "" ] && echo -ne "--dependency $([ "$2" != "" ] && printf "aftercorr" || printf "afterok"):${1}"
	[ "${@}a" != "a" ] && echo -ne "--dependency aftercorr:${@}"
}
export -f depCheckArrack

###################
# Gets job category data
###################
function dispatch {
	jobWait
	JOBCAT=${1}
	printf "%s--account %s --mail-user %s --mail-type %s --time %.0f --mem-per-cpu %d --cpus-per-task %d --profile=task" \
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
echo -e "$HEADER:\tWalltime: $(printHMS $(printSeconds ${SB[$HEADER,MWT]}))"
	echo -e "\tCores:    ${SB[$HEADER,CPT]}"
	echo -e "\tMemory:   $((${SB[$HEADER,MPC]} * ${SB[$HEADER,CPT]}))"
}
export -f jobStats

##################
# Delays job submission based on project's previous submitted job.
#
# Limits submission rate to 1 job every MAX_JOB_RATE seconds.
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
export -f jobWait

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

#######################
# Condenses list of numbers to ranges
#
# 1,2,3,4,7,8,9,12,13,14 -> 1-3,4,7-9,12-14
#######################
function condenseList {
	echo "${@}," | \
		sed "s/,/\n/g" | \
		while read num
		do
		if [[ -z $first ]]
		then
			first=$num
			last=$num
			continue
		fi
		if [[ num -ne $((last + 1)) ]]
		then
			if [[ first -eq last ]]
			then
				echo $first
			else
				echo $first-$last
			fi
			first=$num
			last=$num
		else
			: $((last++))
		fi
	done | paste -sd ","
}
export -f condenseList

#####################
# Expand comma separated list of ranges to individual elements
#
# 1,3-5,8,10-12 -> 1,3,4,5,8,10,11,12
#####################
function expandList {
	for f in ${1//,/ }; do
		if [[ $f =~ - ]]; then
			a+=( $(seq ${f%-*} 1 ${f#*-}) )
		else
			a+=( $f )
		fi  
	done
	
	a=${a[*]}
	a=${a// /,}
	
	echo $a
}
export -f expandList

###################
# Returns full path to specified file.
###################
function realpath {
	echo $(cd $(dirname $1); pwd)/$(basename $1);
}
export -f realpath

function fileSize {
	sizeString=" kMGTEPYZ"
	sizeBlock=0
	readSize=$(ls -la $1 | awk '{print $5}')
	while [ $(echo "$readSize / 1024 > 0" | bc) -eq 1 ]; do
		#printf "%-12s %.0f%s\n" "Read size" $readSize $(echo ${sizeString:${sizeBlock}:1}Bytes | sed -e 's/ //g')
		readSize=$(echo "$readSize / 1024" | bc -l)
		sizeBlock=$((sizeBlock+1))
	done
	echo $(printf "%.0f" $readSize)${sizeString:${sizeBlock}:1}B | sed -e 's/ //g'
}
export -f fileSize
