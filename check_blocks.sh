#!/bin/bash

source /projects/uoo00032/Resources/bin/baserefs.sh

   SAMPLE=${1}
READGROUP=${2}
  READNUM=${3}
    BLOCK=${4}
     NEXT=${5}

curRead1File=blocks/${READGROUP}_R1_${BLOCK}.fastq.gz
curRead2File=blocks/${READGROUP}_R2_${BLOCK}.fastq.gz

IDN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')
DNA=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $2}')
LIB=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $3}')
RUN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $4}')
PLATFORM=Genomic

function depCheck {
	if [ "$1" != "" ]
	then
		echo -ne "--dependency afterok:${1}"
	fi
}

function depDone {
	if [ "$1" != "" ]
	then
		echo -ne "$1"
	else
		echo -ne "Completed!"
	fi
}

function spoolAlign {
	alignOutput=aligned/${BLOCK}.bam
	sortOutput=sorted/${BLOCK}.bam
	
	mkdir -p $(dirname ${alignOutput})
	mkdir -p $(dirname ${sortOutput})
	
	printf "SA: Alignment "
	
	if [ ! -e ${alignOutput}.done ]; then
		DEP_PA=$(sbatch -J PA_${SAMPLE}_${BLOCK} ${SLSBIN}/align.sl ${SAMPLE} ${READGROUP} ${BLOCK} ${alignOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_PA}
		fi
	else
		printf "done "
	fi
	
	printf "%s " "-> SortSAM"
	
	if [ ! -e ${sortOutput}.done ]; then
		DEP_SS=$(sbatch -J SS_${SAMPLE}_${BLOCK} $(depCheck ${DEP_PA}) ${SLSBIN}/sortsam.sl ${SAMPLE} ${alignOutput} ${sortOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_SS" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_SS}
		fi
	else
		printf "done "
	fi
	
	splitArray=""
	
	# Loop though number of contigs in reference sequence.
	# Build list of non-completed contigs blocks.
	for i in $(seq 1 ${NUMCONTIGS}); do
		splitOutput=split/${CONTIGA[$i]}_${BLOCK}.bam
		mkdir -p $(dirname $splitOutput)
		
		if [ ! -e ${splitOutput}.done ]; then
			# This contig block hasn't been split yet.
			splitArray=$(appendList "$splitArray" $i ",")
		fi
	done
	
	printf "%s " "-> ContigSplit"
	
	if [ "$splitArray" != "" ]; then
		# Split elements defined!
		
		DEP_CS=$(sbatch -J CS_${SAMPLE}_${BLOCK} -a $splitArray $(depCheck ${DEP_SS}) ${SLSBIN}/contigsplit.sl ${SAMPLE} ${sortOutput} ${BLOCK} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CS" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s [%02d]" "$DEP_CS" $(splitByChar "$splitArray" "," | wc -w)
		fi
		mkdir -p mergeDeps
		echo "${DEP_CS}" > mergeDeps/merge_${DEP_CS}.sh
	else
		# No array elements define. They're all done!
		printf "done "
	fi
	
	printf "\n"
	
	if [ "$BLOCK" == "$NEXT" ]; then	# This is the final chunk so we can spool up the next step
		echo "CB: Last split so spool up block merging."
		spoolMerge
	fi
}

function spoolMerge {
	READGROUP=$(cat blocks/R1_ReadGroup.txt)
	
	cd ../
	
	mkdir -p slurm
	
	if [ $(echo ${READGROUP} | wc -w) -gt 1 ]
	then
		echo "CB: Too many read-groups. No spooling further pipeline."
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
	numContigs=$(echo $CONTIGS | wc -w)
	echo "CB: Merge Dependencies: ${numDeps}/${numContigs}:[${mergeDeps}]"
	
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
		contig=${CONTIGA[$i]}	# Does bash do array lookups every time too?
		
		printf "%2d:%-10s " $i "$contig"
		
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
			printf "MC "
		fi
		
		if [ ! -e ${markOutput}.done ]; then
			markArray=$(appendList "$markArray"  $i ",")
			printf "MD "
		fi
		
		if [ ! -e ${baseOutput}.done ]; then
			baseArray=$(appendList "$baseArray"  $i ",")
			printf "BR "
		fi
		
		if [ ! -e ${printOutput}.done ]; then
			printArray=$(appendList "$printArray" $i ",")
			printf "PR "
		fi
		
		if [ ! -e ${depthOutput}.done ]; then
			depthArray=$(appendList "$depthArray" $i ",")
			printf "DC "
		fi
		
		if [ "$contig" != "X" ] && [ "$contig" != "Y" ] && [ "$contig" != "MT" ]; then	#Skip sex and mitochondrial chromosomes
			if [ ! -e ${haploOutput}.done ]; then
				haploArray=$(appendList "$haploArray" $i ",")
				printf "HC "
			fi
		fi
		printf "\n"
	done
	
	mergeReadCount=$(echo ${catInputs} | wc -w)
	
	printf "SM: Merge  " 
	
	if [ "$mergeArray" != "" ]; then
		DEP_CM=$(sbatch -J MC_${IDN} -a $mergeArray $(depCheck ${mergeDeps}) ${SLSBIN}/mergecontigs.sl | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CM" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s [%02d] " "$DEP_CM" $(splitByChar "$mergeArray" "," | wc -w)
		fi
	else
		printf "done "
	fi
	
	printf "%s " "-> MarkDup"
	
	if [ "$markArray" != "" ]; then
		DEP_MD=$(sbatch -J MD_${IDN} -a $markArray --begin=now+1hour ${SLSBIN}/markduplicates.sl | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_MD" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s [%02d] " "$DEP_MD" $(splitByChar "$markArray" "," | wc -w)
		fi
		
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$markArray" "$DEP_MD" "$mergeArray" "$DEP_CM"
		scontrol update JobId=${DEP_MD} StartTime=now+0
	else
		printf "done "
	fi
	
	printf "%s " "-> BaseRecal"
	
	if [ "$baseArray" != "" ]; then
		DEP_BR=$(sbatch -J BR_${IDN} -a $baseArray --begin=now+1hour ${SLSBIN}/baserecalibrator.sl | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_BR" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s [%02d] " "$DEP_BR" $(splitByChar "$baseArray" "," | wc -w)
		fi
		
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$baseArray" "$DEP_BR" "$markArray" "$DEP_MD"
		scontrol update JobId=${DEP_BR} StartTime=now+0
	else
		printf "done "
	fi
	
	printf "%s " "-> PrintReads"
	
	if [ "$printArray" != "" ]; then
		DEP_PR=$(sbatch -J PR_${IDN} -a $printArray --begin=now+1hour ${SLSBIN}/printreads.sl | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_PR" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s [%02d] " "$DEP_PR" $(splitByChar "$printArray" "," | wc -w)
		fi
		
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$printArray" "$DEP_PR" "$baseArray" "$DEP_BR"
		scontrol update JobId=${DEP_PR} StartTime=now+0
	else
		printf "done "
	fi
	
	printf "%s " "-> Depth"
	
	if [ "$depthArray" != "" ]; then
		DEP_DC=$(sbatch -J DC_${IDN} -a $depthArray --begin=now+1hour ${SLSBIN}/depthofcoverage.sl ${PLATFORM} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_DC" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s [%02d] " "$DEP_DC" $(splitByChar "$depthArray" "," | wc -w)
		fi
		
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$depthArray" "$DEP_DC" "$printArray" "$DEP_PR"
		scontrol update JobId=${DEP_DC} StartTime=now+0
	else
		printf "done "
	fi
	
	printf "%s " "-> Coverage"
	
	if [ ! -e coverage.sh.done ]; then
		DEP_GD=$(sbatch -J GD_${IDN} $(depCheck ${DEP_DC}) ${SLSBIN}/coverage.sl ${IDN} ${PLATFORM} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_GD" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "$s " "${DEP_GD}"
		fi
	else
		printf "done -> Save Coverage "
		
		# Check if coverage map has been uploaded!
		if [ ! -e ${IDN}.coverage.sh.transfer.done ]; then
			DEP_TGD=$(sbatch -J TGD_${IDN} $(depCheck ${DEP_GD}) ${SLSBIN}/transfer.sl ${IDN} coverage.sh ${IDN}.coverage.sh | awk '{print $4}')
			if [ $? -ne 0 ] || [ "$DEP_TGD" == "" ]; then
				printf "FAILED!\n"
				exit 1
			else
				printf "$s " "${DEP_TGD}"
			fi
		else
			printf "done "
		fi
	fi
	
	printf "%s " "-> Haplo"
	
	if [ "$haploArray" != "" ]; then
		DEP_HC=$(sbatch -J HC_${IDN} -a $haploArray ${SLSBIN}/haplotypecaller.sl | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HC" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s [%02d] " "$DEP_HC" $(splitByChar "$haploArray" "," | wc -w)
		fi
	else
		printf "done "
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
	
	printf "%s " "-> Gender-Haplo"
	
	if [ ! -e ${haploXPar1Output}.done ]; then
		DEP_HCXPAR1=$(sbatch -J HC_${IDN}_${XPAR1} $(depCheck ${DEP_GD}) ${SLSBIN}/haplotypecaller.sl  "${XPAR1}" | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HCXPAR1" == "" ]; then
			printf "%s FAILED!\n" "$XPAR1"
			exit 1
		else
			printf "%s %s " "$XPAR1" "$DEP_HCXPAR1"
		fi
	else
		printf "%s done " "$XPAR1"
	fi
	
	if [ ! -e ${haploTRUEXOutput}.done ]; then
		DEP_HCTRUEX=$(sbatch -J HC_${IDN}_${TRUEX} $(depCheck ${DEP_GD}) ${SLSBIN}/haplotypecaller.sl "${TRUEX}" | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HCTRUEX" == "" ]; then
			printf "%s FAILED!\n" "$TRUEX"
			exit 1
		else
			printf "%s %s " "$TRUEX" "$DEP_HCTRUEX"
		fi
	else
		printf "%s done " "$TRUEX"
	fi
	
	if [ ! -e ${haploXPar2Output}.done ]; then
		DEP_HCXPAR2=$(sbatch -J HC_${IDN}_${XPAR2} $(depCheck ${DEP_GD}) ${SLSBIN}/haplotypecaller.sl "${XPAR2}" | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HCXPAR2" == "" ]; then
			printf "%s FAILED!\n" "$XPAR2"
			exit 1
		else
			printf "%s %s " "$XPAR2" "$DEP_HCXPAR2"
		fi
	else
		printf "%s done " "$XPAR2"
	fi
	
	if [ ! -e ${haploYOutput}.done ]; then
		DEP_HCY=$(sbatch -J HC_${IDN}_Y $(depCheck ${DEP_GD}) ${SLSBIN}/haplotypecaller.sl "Y" | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HCY" == "" ]; then
			printf "Y FAILED!\n"
			exit 1
		else
			printf "Y %s " "$DEP_HCY"
		fi
	else
		printf "Y done"
	fi
	
	# If something went wrong we want to pick up where we left off without stray : characters.
	CatVarDeps="${DEP_HC} ${DEP_HCXPAR1} ${DEP_HCTRUEX} ${DEP_HCXPAR2} ${DEP_HCY}"
	CarVarDeps2=""
	for dep in ${CatVarDeps}; do
		CarVarDeps2=$(appendList "$CarVarDeps2" "${dep}" ":")
	done
	CatVarDeps="${CarVarDeps2}"
	
	CatVarInputs=""
	for contig in ${CONTIGS}; do
		if [ "$contig" == "MT" ]; then
			continue
		elif [ "$contig" == "X" ]; then
			thisInput="haplo/${XPAR1}.g.vcf.gz haplo/${TRUEX}.g.vcf.gz haplo/${XPAR2}.g.vcf.gz"
		else
			thisInput="haplo/${contig}.g.vcf.gz"
		fi
		
		CatVarInputs=$(appendList "${CatVarInputs}" "${thisInput}")
	done
	
	mergePrintOutput=${IDN}.bam
	catvarOutput=${IDN}.g.vcf.gz
	
	printf "\nSM: CatReads "
	
	# Merge print-read bams.
	if [ ! -e ${mergePrintOutput}.done ]; then
		DEP_CR=$(sbatch -J CR_${IDN} $(depCheck ${DEP_PR}) ${SLSBIN}/catreads.sl "${catInputs}" ${mergePrintOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CR" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "$d " ${DEP_CR}
		fi
	else
		printf "done %s " "-> Save Reads"
		
		if [ ! -e ${mergePrintOutput}.transfer.done ]; then
			DEP_TCR=$(sbatch -J TR_${IDN} $(depCheck ${DEP_CR}) ${SLSBIN}/transfer.sl ${IDN} ${mergePrintOutput} | awk '{print $4}')
			if [ $? -ne 0 ] || [ "$DEP_TCR" == "" ]; then
				printf "FAILED!\n"
				exit 1
			else
				printf "$d " ${DEP_TCR}
			fi
		else
			printf "done "
		fi
		
		printf "%s " "-> Reads Index"
		
		if [ ! -e ${mergePrintOutput%.bam}.bai.done ]; then
			DEP_RI=$(sbatch -J RI_${IDN} $(depCheck ${DEP_CR}) ${SLSBIN}/catreadsindex.sl ${mergePrintOutput} ${mergePrintOutput%.bam}.bai | awk '{print $4}')
			if [ $? -ne 0 ] || [ "$DEP_RI" == "" ]; then
				printf "FAILED!\n"
				exit 1
			else
				printf "$d " ${DEP_RI}
			fi
		else
			printf "done %s " "-> Save Reads Index"
			
			if [ ! -e ${mergePrintOutput%.bam}.bai.transfer.done ]; then
				DEP_TRI=$(sbatch -J TRI_${IDN} $(depCheck ${DEP_RI}) ${SLSBIN}/transfer.sl ${IDN} ${mergePrintOutput%.bam}.bai | awk '{print $4}')
				if [ $? -ne 0 ] || [ "$DEP_TRI" == "" ]; then
					printf "FAILED!\n"
					exit 1
				else
					printf "%d " ${DEP_TRI}
				fi
			else
				printf "done "
			fi
		fi
	fi
	
#	printf "\nSM: Fingerprint "
#	
#	if [ ! -e ${IDN}.fingerprint.vcf.gz.done ]; then
#		DEP_FP=$(sbatch -J FP_${IDN} $(depCheck ${DEP_CR}) ${SLSBIN}/fingerprint.sl ${IDN} ${mergePrintOutput} ${IDN}.fingerprint.vcf.gz | awk '{print $4}')
#		if [ $? -ne 0 ] || [ "$DEP_FP" == "" ]; then
#			printf "%s FAILED!\n" "$DEP_FP"
#			exit 1
#		else
#			printf "%s " "$DEP_FP"
#		fi
#	else
#		printf "done "
#	fi
	
	printf "\nSM: CatVariants "
	
	if [ ! -e ${catvarOutput}.done ]; then
		DEP_CV=$(sbatch -J CV_${IDN} $(depCheck ${CatVarDeps}) ${SLSBIN}/catvar.sl "${CatVarInputs}" ${catvarOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CV" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d" ${DEP_CV}
		fi
	else
		printf "done %s " "-> Save Cat Variants" 
		
		if [ ! -e ${catvarOutput}.transfer.done ]; then
			DEP_TV=$(sbatch -J TV_${IDN} $(depCheck ${DEP_CV}) ${SLSBIN}/transfer.sl ${IDN} ${catvarOutput} | awk '{print $4}')
			if [ $? -ne 0 ] || [ "$DEP_TV" == "" ]; then
				printf "FAILED!\n"
				exit 1
			else
				printf "%d " ${DEP_TV}
			fi
		else
			printf "done "
		fi
		
		printf "& Save Variants Index "
		
		if [ ! -e ${catvarOutput}.tbi.transfer.done ]; then
			DEP_TVI=$(sbatch -J TVI_${IDN} $(depCheck ${DEP_CV}) ${SLSBIN}/transfer.sl ${IDN} ${catvarOutput}.tbi | awk '{print $4}')
			if [ $? -ne 0 ] || [ "$DEP_TVI" == "" ]; then
				printf "FAILED!\n"
				exit 1
			else
				printf "%d " ${DEP_TVI}
			fi
		else
			printf "done"
		fi
		
	fi
	
	if [ "$DEP_CV" != "" ] && [ "$DEP_CR" != "" ]; then
		saveDeps="${DEP_CV}:${DEP_CR}"
	elif [ "$DEP_CR" != "" ]; then
		saveDeps="${DEP_CR}"
	elif [ "$DEP_CV" != "" ]; then
		saveDeps="${DEP_CV}"
	fi
	
	printf "SM: Save metrics "
	
	if [ ! -e ${IDN}.metrics.txt.transfer.done ]; then
		DEP_TMet=$(sbatch -J TMet_${IDN} $(depCheck ${saveDeps}) ${SLSBIN}/transfer.sl ${IDN} metrics.txt ${IDN}.metrics.txt | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_TMet" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_TMet}
		fi
	else
		printf "done "
	fi
	
	printf "& Save commands "
	
	if [ ! -e ${IDN}.commands.txt.transfer.done ]; then
		DEP_TCMD=$(sbatch -J TCMD_${IDN} $(depCheck ${saveDeps}) ${SLSBIN}/transfer.sl ${IDN} commands.txt ${IDN}.commands.txt | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_TCMD" == "" ]; then
			printf "FAILED!"
			exit 1
		else
			printf "%d" ${DEP_TCMD}
		fi
	else
		printf "done "
	fi
	
	printf "\n"
}

if [ "$READNUM" == "R1" ]; then	# BLOCK Read1 has completed!
	if [ -e ${curRead2File}.done ]; then	# NEXT Read2 exists so BLOCK Read2 has completed!
		echo "CB: Both R1 and R2 ${BLOCK} completed!" #| tee -a check_${BLOCK}.txt
		spoolAlign
	else	# NEXT Read2 doesn't exist yet so BLOCK Read2 hasn't completed!
		echo "CB: R1 ${BLOCK} completed. Waiting for R2!" #| tee -a check_${BLOCK}.txt
	fi
elif [ "$READNUM" == "R2" ]; then	# BLOCK Read2 has completed!
	if [ -e ${curRead1File}.done ]; then		# NEXT Read1 exists so BLOCK Read1 has completed!
		echo "CB: Both R1 and R2 ${BLOCK} completed!" #| tee -a check_${BLOCK}.txt
		spoolAlign
	else	# NEXT Read1 doesn't exist yet so BLOCK Read1 hasn't completed!
		echo "CB: R2 ${BLOCK} complete. Waiting for R1!" #| tee -a check_${BLOCK}.txt
	fi
fi
