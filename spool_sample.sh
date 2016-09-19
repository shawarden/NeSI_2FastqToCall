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

printf "%-22s%s\n" "SampleID" "${SAMPLE}"
printf "%-22s%s\n" "Read 1" "${READ1}"
printf "%-22s%s\n" "Read 2" "${READ2}"
printf "%-22s%s\n" "PLATFORM" "${PLATFORM}"

IDN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')
DNA=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $2}')
LIB=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $3}')
RUN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $4}')

# Platform setup.
platformBED=${PLATFORMS}/${PLATFORM}.bed
  genderBED=${PLATFORMS}/$([ "${PLATFORM}" == "Genomic" ] && echo "AV5" || echo "${PLATFORM}" ).bed
  genderSRC=${genderBED%.bed}.sh
  
   CUR_PATH=$(pwd)
if [ "$LOCATION" == "scratch" ]; then
	WORK_PATH=/scratch/jobs/$USER
else
	WORK_PATH=/projects/uoo00032/Alignments
fi


SAMPLE_PATH=${WORK_PATH}/${IDN}/${DNA}_${LIB}_${RUN}

if ! mkdir -p ${SAMPLE_PATH}/slurm; then
	echo "Error creating output folder!"
	exit 1
fi

date '+%Y%m%d_%H%M%S' > ${WORK_PATH}/${IDN}/starttime.txt

cd ${SAMPLE_PATH}

alignArray=""

for i in $(seq 0 $FASTQ_MAXJOBZ); do
	curBlock=$(printf "%0${FASTQ_MAXZPAD}d" $i)
	alignOutput=aligned/$curBlock.bam
	sortOutput=sorted/$curBlock.bam
	
	if [ ! -e ${alignOutput}.done ]; then
		# This contig block hasn't been split yet.
		alignArray=$(appendList "$alignArray" $i ",")
	fi
	
	if [ ! -e ${sortOutput}.done ]; then
		# This contig block hasn't been split yet.
		sortArray=$(appendList "$sortArray" $i ",")
	fi
done

# Dispatch alignemnt array.
# Alignemnt array doesn't have a dependency yet since ReadSplit needs to know what the aligner's JobID is to update it. Delay the start.
printf "%-22s" "Align->Sort->Split"
DEP_BA=$(sbatch $(dispatch "BA") -J BA_${SAMPLE} --begin=now+72hour --array=$alignArray ${SLSBIN}/blockalign.sl ${SAMPLE} | awk '{print $4}')
if [ $? -ne 0 ] || [ "$DEP_BA" == "" ]; then
	printf "FAILED!\n"
	exit 1
else
	printf "%sx%-4d\n" "${DEP_BA}" $(splitByChar "$alignArray" "," | wc -w)
fi

# Generate entire job tree now, waiting for AlignBlock to start.

# Things HaplotypeCaller is dependent on.
# This only applied to Chr X and Y and their PAR1/2 regions.
genderDeps=""

# Things CatVariants is dependant on.
# All HaplotypeCaller jobs.
CatVarDeps=""
CatVarInputs=""

mergeArray=""
markArray=""
baseArray=""
printArray=""
depthArray=""
haploArray=""

# Loop though number of contigs in reference sequence.
# Build list of incomplete merged contigs.
for i in $(seq 1 ${NUMCONTIGS}); do
	# Build input/output file names
	contig=${CONTIGARRAY[$i]}	# Does bash do array lookups every time too?
	
	#printf "%2d:%-10s " $i "$contig"
	
	# Gather merge inputs
	mergeOutput=merged/${contig}.bam
	 markOutput=markdup/${contig}.bam
	 baseOutput=baserecal/${contig}.firstpass
	printOutput=printreads/${contig}.bam
	  catInputs=$(appendList "$catInputs" "${printOutput}" " ")
	depthOutput=depth/${contig} #.sample_summary
	haploOutput=haplo/${contig}.g.vcf.gz
	
	mkdir -p $(dirname $mergeOutput)
	mkdir -p $(dirname $markOutput)
	mkdir -p $(dirname $baseOutput)
	mkdir -p $(dirname $printOutput)
	mkdir -p $(dirname $depthOutput)
	mkdir -p $(dirname $haploOutput)
	
	if [ ! -e ${mergeOutput}.done ]; then
		mergeArray=$(appendList "$mergeArray" $i ",")
		#printf "MC "
	fi
	
	if [ ! -e ${markOutput}.done ]; then
		markArray=$(appendList "$markArray"  $i ",")
		#printf "MD "
	fi
	
	if [ ! -e ${baseOutput}.done ]; then
		baseArray=$(appendList "$baseArray"  $i ",")
		#printf "BR "
	fi
	
	if [ ! -e ${printOutput}.done ]; then
		printArray=$(appendList "$printArray" $i ",")
		#printf "PR "
	fi
	
	if [ ! -e ${depthOutput}.done ]; then
		depthArray=$(appendList "$depthArray" $i ",")
		#printf "DC "
	fi
	
	if [ "$contig" != "X" ] && [ "$contig" != "Y" ] && [ "$contig" != "MT" ]; then	#Skip sex and mitochondrial chromosomes
		if [ ! -e ${haploOutput}.done ]; then
			haploArray=$(appendList "$haploArray" $i ",")
			#printf "HC "
		fi
	fi
	#printf "\n"
done

mergeReadCount=$(echo ${catInputs} | wc -w)

printf "%-22s" "Merge Contig" 

if [ "$mergeArray" != "" ]; then
	DEP_CM=$(sbatch $(dispatch "MC") -J MC_${IDN} --array $mergeArray $(depCheck ${DEP_BA}) ${SLSBIN}/mergecontigs.sl | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_CM" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%sx%-2d [%s]\n" "$DEP_CM" $(splitByChar "$mergeArray" "," | wc -w) $(condenseList "$mergeArray")
	fi
else
	printf "done\n"
fi

printf "%-22s" "Mark Duplicates"

if [ "$markArray" != "" ]; then
	DEP_MD=$(sbatch $(dispatch "MD") -J MD_${IDN} --array $markArray $(depCheck ${DEP_CM}) ${SLSBIN}/markduplicates.sl | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_MD" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$markArray" "$DEP_MD" "$mergeArray" "$DEP_CM"
		#scontrol update JobId=${DEP_MD} StartTime=now+0
		printf "%sx%-2d [%s]\n" "$DEP_MD" $(splitByChar "$markArray" "," | wc -w) $(condenseList "$markArray")
	fi
else
	printf "done\n"
fi

printf "%-22s" "Base Recalibration"

if [ "$baseArray" != "" ]; then
	DEP_BR=$(sbatch $(dispatch "BR") -J BR_${IDN} --array $baseArray $(depCheck ${DEP_MD}) ${SLSBIN}/baserecalibrator.sl | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_BR" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$baseArray" "$DEP_BR" "$markArray" "$DEP_MD"
		#scontrol update JobId=${DEP_BR} StartTime=now+0
		printf "%sx%-2d [%s]\n" "$DEP_BR" $(splitByChar "$baseArray" "," | wc -w) $(condenseList "$baseArray")
	fi
else
	printf "done\n"
fi

printf "%-22s" "Print Reads"

if [ "$printArray" != "" ]; then
	DEP_PR=$(sbatch $(dispatch "PR") -J PR_${IDN} --array $printArray $(depCheck ${DEP_BR}) ${SLSBIN}/printreads.sl | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_PR" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$printArray" "$DEP_PR" "$baseArray" "$DEP_BR"
		#scontrol update JobId=${DEP_PR} StartTime=now+0
		printf "%sx%-2d [%s]\n" "$DEP_PR" $(splitByChar "$printArray" "," | wc -w) $(condenseList "$printArray")
	fi
else
	printf "done\n"
fi

printf "%-22s" "Depth of Coverage"

if [ "$depthArray" != "" ]; then
	DEP_DC=$(sbatch $(dispatch "DC") -J DC_${IDN} --array $depthArray $(depCheck ${DEP_PR}) ${SLSBIN}/depthofcoverage.sl ${PLATFORM} | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_DC" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$depthArray" "$DEP_DC" "$printArray" "$DEP_PR"
		#scontrol update JobId=${DEP_DC} StartTime=now+0
		printf "%sx%-2d [%s]\n" "$DEP_DC" $(splitByChar "$depthArray" "," | wc -w) $(condenseList "$depthArray")
	fi
else
	printf "done\n"
fi

printf "%-22s" "Gender Determination"

if [ ! -e coverage.sh.done ]; then
	DEP_GD=$(sbatch $(dispatch "GD") -J GD_${IDN} $(depCheck ${DEP_DC}) ${SLSBIN}/coverage.sl ${IDN} ${PLATFORM} | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_GD" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "${DEP_GD}"
	fi
else
	printf "done\n"
fi

printf "%-22s" "HaplotypeCaller"

if [ "$haploArray" != "" ]; then
	DEP_HC=$(sbatch $(dispatch "HC") -J HC_${IDN} --array $haploArray $(depCheck ${DEP_PR}) ${SLSBIN}/haplotypecaller.sl | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_HC" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$haploArray" "$DEP_HC" "$depthArray" "$DEP_PR"
		#scontrol update JobId=${DEP_HC} StartTime=now+0
		printf "%sx%-2d [%s]\n" "$DEP_HC" $(splitByChar "$haploArray" "," | wc -w) $(condenseList "$haploArray")
	fi
else
	printf "done\n"
fi

CatVarDeps=$(appendList "$CatVarDeps" "${DEP_HC}")

haploXInput=printreads/X.bam
haploYInput=printreads/Y.bam

haploXPar1Output=haplo/${XPAR1}.g.vcf.gz
haploTRUEXOutput=haplo/${TRUEX}.g.vcf.gz
haploXPar2Output=haplo/${XPAR2}.g.vcf.gz
	haploYOutput=haplo/Y.g.vcf.gz
mkdir -p $(dirname ${haploXPar1Output})
mkdir -p $(dirname ${haploTRUEXOutput})
mkdir -p $(dirname ${haploXPar2Output}) 
mkdir -p $(dirname ${haploYOutput})

printf "%-22s" "HaplotypeCaller XPAR1"

if [ ! -e ${haploXPar1Output}.done ]; then
	DEP_HCXPAR1=$(sbatch $(dispatch "HC") -J HC_${IDN}_XPAR1 --array=23 $(depCheck ${DEP_GD}) ${SLSBIN}/haplotypecaller.sl  "${XPAR1}" | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_HCXPAR1" == "" ]; then
		printf "FAILED!\n" 
		exit 1
	else
		printf "%s\n" "$DEP_HCXPAR1"
	fi
else
	printf "done\n"
fi

printf "%-22s" "HaplotypeCaller TRUEX"

if [ ! -e ${haploTRUEXOutput}.done ]; then
	DEP_HCTRUEX=$(sbatch $(dispatch "HC") -J HC_${IDN}_TRUEX --array=23 $(depCheck ${DEP_GD}) ${SLSBIN}/haplotypecaller.sl "${TRUEX}" | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_HCTRUEX" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "$DEP_HCTRUEX"
	fi
else
	printf "done\n"
fi

printf "%-22s" "HaplotypeCaller XPAR2"

if [ ! -e ${haploXPar2Output}.done ]; then
	DEP_HCXPAR2=$(sbatch $(dispatch "HC") -J HC_${IDN}_XPAR2 --array=23 $(depCheck ${DEP_GD}) ${SLSBIN}/haplotypecaller.sl "${XPAR2}" | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_HCXPAR2" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "$DEP_HCXPAR2"
	fi
else
	printf "done\n"
fi

printf "%-22s" "HaplotypeCaller Y"

if [ ! -e ${haploYOutput}.done ]; then
	DEP_HCY=$(sbatch $(dispatch "HC") -J HC_${IDN}_Y --array=24 $(depCheck ${DEP_GD}) ${SLSBIN}/haplotypecaller.sl "Y" | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_HCY" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "$DEP_HCY"
	fi
else
	printf "done\n"
fi

# If something went wrong we want to pick up where we left off without stray : characters.
CatVarDeps="${DEP_HC} ${DEP_HCXPAR1} ${DEP_HCTRUEX} ${DEP_HCXPAR2} ${DEP_HCY}"
CarVarDeps2=""
for dep in ${CatVarDeps}; do
	CarVarDeps2=$(appendList "$CarVarDeps2" "${dep}" ":")
done
CatVarDeps="${CarVarDeps2}"

CatVarInputs=""
for contig in ${CONTIGARRAY[@]}; do
	if [ "$contig" == "MT" ]; then
		continue
	elif [ "$contig" == "X" ]; then
		thisInput="haplo/${XPAR1}.g.vcf.gz haplo/${TRUEX}.g.vcf.gz haplo/${XPAR2}.g.vcf.gz"
	else
		thisInput="haplo/${contig}.g.vcf.gz"
	fi
	
	CatVarInputs=$(appendList "${CatVarInputs}" "${thisInput}")
done

catReadsOutput=${IDN}.bam
catVarOutput=${IDN}.g.vcf.gz

printf "%-22s" "CatReads"

# Merge print-read bams.
if [ ! -e ${catReadsOutput}.done ]; then
	DEP_CR=$(sbatch $(dispatch "CR") -J CR_${IDN} $(depCheck ${DEP_PR}) ${SLSBIN}/catreads.sl "${catInputs}" ${catReadsOutput} | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_CR" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "${DEP_CR}"
	fi
else
	printf "done\n"
	
	printf "%-22s" "Reads Index"
	
	if [ ! -e ${catReadsOutput%.bam}.bai.done ]; then
		DEP_RI=$(sbatch $(dispatch "RI") -J RI_${IDN} $(depCheck ${DEP_CR}) ${SLSBIN}/catreadsindex.sl ${catReadsOutput} ${catReadsOutput%.bam}.bai | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_RI" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s\n" "${DEP_RI}"
		fi
	else
		printf "done\n"
	fi
fi

printf "%-22s" "CatVariants"

if [ ! -e ${catVarOutput}.done ]; then
	DEP_CV=$(sbatch $(dispatch "CV") -J CV_${IDN} $(depCheck ${CatVarDeps}) ${SLSBIN}/catvar.sl "${CatVarInputs}" ${catVarOutput} | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_CV" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "${DEP_CV}"
	fi
else
	printf "done\n" 
fi

if [ "$DEP_CV" != "" ] && [ "$DEP_CR" != "" ]; then
	saveDeps="${DEP_CV}:${DEP_CR}"
elif [ "$DEP_CR" != "" ]; then
	saveDeps="${DEP_CR}"
elif [ "$DEP_CV" != "" ]; then
	saveDeps="${DEP_CV}"
fi

printf "%-22s" "Save metrics"

if [ ! -e ${IDN}.metrics.txt.transfer.done ]; then
	DEP_TMet=$(sbatch $(dispatch "TF") -J TMet_${IDN} $(depCheck ${saveDeps}) ${SLSBIN}/transfer.sl ${IDN} metrics.txt ${IDN}.metrics.txt | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_TMet" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "${DEP_TMet}"
	fi
else
	printf "done\n"
fi

printf "%-22s" "Save commands"

if [ ! -e ${IDN}.commands.txt.transfer.done ]; then
	DEP_TCMD=$(sbatch $(dispatch "TF") -J TCMD_${IDN} $(depCheck ${saveDeps}) ${SLSBIN}/transfer.sl ${IDN} commands.txt ${IDN}.commands.txt | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_TCMD" == "" ]; then
		printf "FAILED!"
		exit 1
	else
		printf "%s\n" "${DEP_TCMD}"
	fi
else
	printf "done\n"
fi

printf "%-22s" "ReadSplitter"

# Split read 1 and 2 into chunks
splitReadArray=""

if [ ! -e ${SAMPLE}_R1_split.done ]; then
	splitReadArray=$(appendList "$splitReadArray" 1 ",")
fi

if [ ! -e ${SAMPLE}_R2_split.done ]; then
	splitReadArray=$(appendList "$splitReadArray" 2 ",")
fi

# Read 1 split isn't complete, run it now.

# Get a fancy size
sizeString=" kMGTEPYZ"
sizeBlock=0
readSize=$(($(ls -la ${READ1} | awk '{print $5}') + $(ls -la ${READ2} | awk '{print $5}')))
while [ $(echo "$readSize / 1024 > 0" | bc) -eq 1 ]; do
	#printf "%-12s %.0f%s\n" "Read size" $readSize $(echo ${sizeString:${sizeBlock}:1}Bytes | sed -e 's/ //g')
	readSize=$(echo "$readSize / 1024" | bc -l)
	sizeBlock=$((sizeBlock+1))
done

readSize=$(echo $(printf "%.0f" $readSize)${sizeString:${sizeBlock}:1}B | sed -e 's/ //g')
#printf "%-12s $s\n" "Read size" $readSize

if [ "$splitReadArray" != "" ]; then
	DEP_SR=$(sbatch $(dispatch "RS") -J RS_${SAMPLE}_${readSize} -a $splitReadArray ${SLSBIN}/readsplit.sl ${SAMPLE} ${READ1} ${READ2} ${DEP_BA} | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_SR" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
	printf "%sx%-1d [%s]\n" "${DEP_SR}" $(splitByChar "$splitReadArray" "," | wc -w) "$splitReadArray"
		
		# Now that we have the ReadSplit jobID, assign the entire Alignment array to it.
		scontrol update JobID=${DEP_BA} StartTime=now Dependency=afterok:${DEP_SR}
	fi
else
	printf "done -> Aligning..."
	# Both reads are done. Build read list and spool alignments.
	READGROUP=$(cat blocks/R1_ReadGroup.txt)
	readBlocks=$(($(find ./blocks -type f -iname "${READGROUP}_R1_*.fastq.gz.done" | wc -l) - 1))
	purgeList="$(($readBlocks+1))-$FASTQ_MAXJOBZ"
	scancel ${DEP_BA}_[${purgeList}] && echo "Purged excess alignment and sort jobs $purgeList"
	for i in $(seq 0 ${readBlocks}); do
		if [ $i -lt $readBlocks ]; then
			nextBlock=$(($i + 1))
		else
			nextBlock=$i
		fi
		${PBIN}/check_blocks.sh ${SAMPLE} ${READGROUP} R1 $i $nextBlock $DEP_BA $DEP_SS
	done
fi


