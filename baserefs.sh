#!/bin/bash

###########
# Base Reference file contains everything I don't want to retype in every script.
#
# include: "source /path/to/baserefs.sh" at the top of pretty much every script to import these values.
###########

# General References.
export      PROJECT=/projects/uoo00032
export    RESOURCES=${PROJECT}/Resources
export    PLATFORMS=${RESOURCES}/Capture_Platforms/GRCh37
export         PBIN=${RESOURCES}/bin
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
export JOB_TEMP_DIR=${TMPDIR}	#system defined temp dir. /tmp/jobs/$SLURM_JOB_USER/$SLURM_JOB_ID

export      BWA=${RESOURCES}/bwa-0.7.15/bwa
export   PICARD=${RESOURCES}/picard-tools-2.5.0/picard.jar
export SAMTOOLS=${RESOURCES}/samtools-1.3.1/samtools
export     GATK=${RESOURCES}/GenomeAnalysisTK-nightly-2016-07-22-g8923599/GenomeAnalysisTK.jar

####################
# Modules versions #
####################

export     MOD_ZLIB="zlib"
export     MOD_JAVA="Java/1.8.0_5"

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

export FASTQ_MAXREAD=10000000	# How many reads per block.
export FASTQ_MAXSCAN=10000		# How many lines to check for best index.
export FASTQ_MAXDIFF=2			# Maximum index variation before new index is created.

###########################
# Java memory calculation #
###########################

# Sometimes we call this file outside of a slurm script so fake it til it's made.
[ -z $SLURM_MEM_PER_CPU ] && SLURM_MEM_PER_CPU=4096
[ -z $SLURM_JOB_CPUS_PER_NODE ] && SLURM_JOB_CPUS_PER_NODE=4

export JAVA_MEM_GB=$((($SLURM_MEM_PER_CPU * $SLURM_JOB_CPUS_PER_NODE)/1024))
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

#######################
# Output basic node information on job failure.
#######################
function scriptFailed {
	HEADER=${1}
	echo ""
	echo "${HEADER}: SControl -----"
	scontrol show job ${SLURM_JOBID}
	echo ""
	echo "${HEADER}: Export -----"
	export
	echo ""
	echo "${HEADER}: Storage -----"
	df -ah
	echo ""
}

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
		"%19s %-50s %-5d %-6d %s\n" \
		"$(date '+%Y-%m-%d %H:%M:%S')" \
		${SLURM_JOB_NAME}$([ "$SLURM_ARRAY_TASK_ID" != "" ] && printf "_%s" "${CONTIGA[${SLURM_ARRAY_TASK_ID}]}") \
		${SLURM_JOB_CPUS_PER_NODE} \
		$((${SLURM_JOB_CPUS_PER_NODE} * ${SLURM_MEM_PER_CPU})) \
		$(printHMS $SECONDS) | \
		tee -a ${BACKDIR}metrics.txt >> ${HOME}/metrics.txt
}

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

# Export all function for use in calling script.
export -f printHMS
export -f scriptFailed
export -f storeMetrics
export -f appendList
export -f splitByChar
export -f tieTaskDeps