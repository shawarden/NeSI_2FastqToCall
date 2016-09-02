#!/bin/bash

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

   SAMPLE=${1}
READGROUP=${2}
  READNUM=${3}
    BLOCK=${4}
     NEXT=${5}
 ALIGNARR=${6}
  SORTARR=${7}
  
ZPADBLOCK=$(printf "%0${FASTQ_MAXZPAD}d" $BLOCK)
#ZPADNEXT=$(printf "%0${FASTQ_MAXZPAD}d" $NEXT)

curRead1File=blocks/${READGROUP}_R1_${ZPADBLOCK}.fastq.gz
curRead2File=blocks/${READGROUP}_R2_${ZPADBLOCK}.fastq.gz

IDN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')
DNA=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $2}')
LIB=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $3}')
RUN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $4}')
PLATFORM=Genomic

function spoolAlign {
	alignOutput=aligned/${ZPADBLOCK}.bam
	sortOutput=sorted/${ZPADBLOCK}.bam
	
	mkdir -p $(dirname ${alignOutput})
	mkdir -p $(dirname ${sortOutput})
	
	printf "SA: Alignment   "
	
	if [ ! -e ${alignOutput}.done ]; then
		#DEP_PA=$(sbatch -J PA_${SAMPLE}_${ZPADBLOCK} ${SLSBIN}/align.sl ${SAMPLE} ${READGROUP} ${ZPADBLOCK} ${alignOutput} | awk '{print $4}')
		# Remove start delay or dependency on this job as the inputs now exist!
		scontrol update JobId=${ALIGNARR}_${BLOCK} StartTime=now Dependency=
		if [ $? -ne 0 ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s\n" ${ALIGNARR}_${BLOCK}
		fi
	else
		printf "done\n"
	fi
	
	printf "SA: SortSAM     "
	
	if [ ! -e ${sortOutput}.done ]; then
		#DEP_SS=$(sbatch -J SS_${SAMPLE}_${ZPADBLOCK} $(depCheck ${DEP_PA}) ${SLSBIN}/sortsam.sl ${SAMPLE} ${alignOutput} ${sortOutput} | awk '{print $4}')
		# Set this element to be dependent on alignments.
		scontrol update JobId=${SORTARR}_${BLOCK} Dependency=${ALIGNARR}_${BLOCK}
		if [ $? -ne 0 ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s\n" ${SORTARR}_${BLOCK}
		fi
	else
		printf "done\n"
	fi
	
	printf "SA: ContigSplit "
	
	splitArray=""
	
	# Loop though number of contigs in reference sequence.
	# Build list of non-completed contigs blocks.
	for i in $(seq 1 ${NUMCONTIGS}); do
		splitOutput=split/${CONTIGARRAY[$i]}_${ZPADBLOCK}.bam
		mkdir -p $(dirname $splitOutput)
		
		if [ ! -e ${splitOutput}.done ]; then
			# This contig block hasn't been split yet.
			splitArray=$(appendList "$splitArray" $i ",")
		fi
	done
	
	if [ "$splitArray" != "" ]; then
		# Split elements defined!
		
		DEP_CS=$(sbatch $(dispatch "CS") -J CS_${SAMPLE}_${BLOCK} --array $splitArray $(depCheck ${SORTARR}_${BLOCK}) ${SLSBIN}/contigsplit.sl ${sortOutput} ${ZPADBLOCK} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CS" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%sx%-2d [%s]\n" "$DEP_CS" $(splitByChar "$splitArray" "," | wc -w) $(condenseList "$splitArray")
		fi
		mkdir -p mergeDeps
		echo "${DEP_CS}" > mergeDeps/merge_${DEP_CS}.sh
	else
		# No array elements define. They're all done!
		printf "done\n"
	fi
	
	if [ "$BLOCK" == "$NEXT" ]; then	# This is the final chunk so we can spool up the next step
		echo "CB: $BLOCK of $NEXT block$([ $NEXT -gt 1 ] && echo -ne "s") completed!"
		purgeList="$((${NEXT}+1))-$FASTQ_MAXJOBZ"
		# Purge extra align and sort array elements.
		scancel ${ALIGNARR}_[${purgeList}] ${SORTARR}_[${purgeList}] && echo "CB: Purged ${ALIGNARR}_[${purgeList}] and ${SORTARR}_[${purgeList}]"
		
		# Spool second segment.
		echo "CB: Spooling merge segement!"
		spoolMerge
	fi
}

function spoolMerge {
	READGROUP=$(cat blocks/R1_ReadGroup.txt)
	
	cd ../
	
	mkdir -p slurm
	
	if [ $(echo ${READGROUP} | wc -w) -gt 1 ]
	then
		echo "CB: Too many read-groups. Not spooling merge and call pipeline."
		exit 1
	fi
	
	# List of printRead jobs that need to complete before picard MergeSamFiles and global Depth of Coverage can occur.
	mergeReadDeps=""
	mergeReadInputs=""
	
	# Things HaplotypeCaller is dependent on.
	# This only applied to Chr X and Y and their PAR1/2 regions.
	genderDeps=""
	
	# Things CatVariants is dependant on.
	# All HaplotypeCaller jobs.
	CatVarDeps=""
	CatVarInputs=""
	
	mergeDeps=""
	mergeDepFiles=$(find . -type f -iwholename "*/mergeDeps/merge_*.sh")
	for file in $mergeDepFiles; do
		mergeDeps=$(appendList "$mergeDeps" "$(cat $file)" ":")
		rm $file
	done
	numDeps=$(echo $mergeDepFiles | wc -w)
	echo "CB: Merge Dependencies: ${numDeps}/${NEXT}:[${mergeDeps}]"
	
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
		
		#printf "SM: %2d:%-10s " $i "$contig"
		
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
	
	printf "SM: %-22s" "Merge Contig" 
	
	if [ "$mergeArray" != "" ]; then
		DEP_CM=$(sbatch $(dispatch "MC") -J MC_${IDN} --array $mergeArray $(depCheck ${mergeDeps}) ${SLSBIN}/mergecontigs.sl | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CM" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%sx%-2d [%s]\n" "$DEP_CM" $(splitByChar "$mergeArray" "," | wc -w) $(condenseList "$mergeArray")
		fi
	else
		printf "done\n"
	fi
	
	printf "SM: %-22s" "Mark Duplicates"
	
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
	
	printf "SM: %-22s" "Base Recalibration"
	
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
	
	printf "SM: %-22s" "Print Reads"
	
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
	
	printf "SM: %-22s" "Depth of Coverage"
	
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
	
	printf "SM: %-22s" "Gender Determination"
	
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
	
	printf "SM: %-22s" "HaplotypeCaller"
	
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
	
	printf "SM: %-22s" "HaplotypeCaller XPAR1"
	
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
	
	printf "SM: %-22s" "HaplotypeCaller TRUEX"
	
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
	
	printf "SM: %-22s" "HaplotypeCaller XPAR2"
	
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
	
	printf "SM: %-22s" "HaplotypeCaller Y"
	
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
	
	printf "SM: %-22s" "CatReads"
	
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
		
		printf "SM: %21s" "Reads Index"
		
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
	
#	printf "SM: %-22s" "Fingerprint"
#	
#	if [ ! -e ${IDN}.fingerprint.vcf.gz.done ]; then
#		DEP_FP=$(sbatch -J FP_${IDN} $(depCheck ${DEP_PR}) ${SLSBIN}/fingerprint.sl ${IDN} fp.haplotyped.vcf.gz | awk '{print $4}')
#		if [ $? -ne 0 ] || [ "$DEP_FP" == "" ]; then
#			printf "%s FAILED!\n" "$DEP_FP"
#			exit 1
#		else
#			printf "%s " "$DEP_FP"
#		fi
#	else
#		printf "done "
#	fi
	
	printf "SM: %-22s" "CatVariants"
	
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
	
	printf "SM: %-22s" "Save metrics"
	
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
	
	printf "SM: %-22s" "Save commands"
	
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
}

if [ "$READNUM" == "R1" ]; then	# BLOCK Read1 has completed!
	if [ -e ${curRead2File}.done ]; then	# NEXT Read2 exists so BLOCK Read2 has completed!
		echo "CB: Both R1 and R2 ${BLOCK} completed!" #| tee -a check_${BLOCK}.txt
		spoolAlign
	else	# NEXT Read2 doesn't exist yet so BLOCK Read2 hasn't completed!
		echo "CB: R1 ${BLOCK} completed but R2 not done yet!" #| tee -a check_${BLOCK}.txt
	fi
elif [ "$READNUM" == "R2" ]; then	# BLOCK Read2 has completed!
	if [ -e ${curRead1File}.done ]; then		# NEXT Read1 exists so BLOCK Read1 has completed!
		echo "CB: Both R1 and R2 ${BLOCK} completed!" #| tee -a check_${BLOCK}.txt
		spoolAlign
	else	# NEXT Read1 doesn't exist yet so BLOCK Read1 hasn't completed!
		echo "CB: R2 ${BLOCK} complete but R1 not done yet!" #| tee -a check_${BLOCK}.txt
	fi
fi
