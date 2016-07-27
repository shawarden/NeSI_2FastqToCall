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

export CONTIGS=$(cat ${REFD} | awk -F'[[:space:]:]' 'NR!=1{print $3}')
export NUMCONTIGS=$(echo "${CONTIGS}" | wc -w)

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

######################
# SBATCH definitions #
######################
declare -A SB

# SB[W,*]: sbatch max wall time in minutes. Keep this at least 1.25x highest run-time.
# SB[C,*]: sbatch cores per task
# SB[M,*]: sbatch memory per task in megabytes

SB[ACCOUNT]=uoo00032
SB[MAIL,USER]=sam.hawarden@otago.ac.nz
SB[MAIL,TYPE]=FAIL

SB[W,READSPLIT]=300
SB[C,READSPLIT]=8
SB[M,READSPLIT]=16384

SB[W,ALIGN]=120
SB[C,ALIGN]=6
SB[M,ALIGN]=16384

SB[W,SORTSAM]=60
SB[C,SORTSAM]=4
SB[M,SORTSAM]=16384

SB[W,CONTIGSPLIT]=60
SB[C,CONTIGSPLIT]=2
SB[M,CONTIGSPLIT]=4096

SB[W,MERGECONTIG]=90
SB[C,MERGECONTIG]=4
SB[M,MERGECONTIG]=16384

SB[W,MARKDUP]=120
SB[C,MARKDUP]=2
SB[M,MARKDUP]=16384

SB[W,BASERECAL]=120
SB[C,BASERECAL]=4
SB[M,BASERECAL]=16384

SB[W,PRINT]=180
SB[C,PRINT]=4
SB[M,PRINT]=16384

SB[W,HAPLO]=180
SB[C,HAPLO]=8
SB[M,HAPLO]=32768

# Contig Wall time multipliers. Do we really need 84 of these? up. Move to another file!
SB[CWM,1]=1
SB[CWM,2]=1
SB[CWM,3]=1
SB[CWM,4]=1
SB[CWM,5]=1
SB[CWM,6]=1

SB_ARGS="--account ${SB[ACCOUNT]} --mail-user=${SB[MAIL,USER]}"

export OPENBLAS_MAIN_FREE=1

####################
# Helper functions #
####################

function printHMS {
	SECS=${1}
	printf "%02d:%02d:%02d" "$(($SECS / 3600))" "$((($SECS % 3600) / 60))" "$(($SECS % 60))"
}

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

function storeMetrics {
	if [ "${1}" != "" ]; then
		BACKDIR="../"
	else
		BACKDIR=""
	fi
	
	printf \
	"%19s %-50s %-5d %-6d %02d:%02d:%02d\n" \
	"$(date '+%Y-%m-%d %H:%M:%S')" \
	"$SLURM_JOB_NAME" \
	"$SLURM_JOB_CPUS_PER_NODE" \
	"$(($SLURM_JOB_CPUS_PER_NODE * $SLURM_MEM_PER_CPU))" \
	"$(($SECONDS / 3600))" \
	"$((($SECONDS % 3600) / 60))" \
	"$(($SECONDS % 60))" | \
	tee -a ${BACKDIR}metrics.txt >> $HOME/metrics.txt
}

export -f printHMS
export -f scriptFailed
export -f storeMetrics