#!/bin/bash

source /projects/uoo00032/Resources/bin/baserefs.sh

   SAMPLE=${1}
READGROUP=${2}
  READNUM=${3}
    BLOCK=${4}
     NEXT=${5}

curRead1File=blocks/${READGROUP}_R1_${BLOCK}.fastq.gz
curRead2File=blocks/${READGROUP}_R2_${BLOCK}.fastq.gz

mkdir -p aligned
mkdir -p sorted

IDN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')
DNA=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $2}')
LIB=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $3}')
RUN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $4}')
PLATFORM=Genomic

# Platform setup.
platformBED=${PLATFORMS}/${PLATFORM}.bed
  genderBED=${PLATFORMS}/$([ "${PLATFORM}" == "Genomic" ] && echo "AV5" || echo "${PLATFORM}" ).bed
  genderSRC=${genderBED%.bed}.sh
  
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
	alignOutput=aligned/${READGROUP}_${BLOCK}.bam
	sortOutput=sorted/${READGROUP}_${BLOCK}.bam
	
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
		splitOutput=split/${CONTIGA[$i]}/split_${BLOCK}.bam
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
			printf "%d [%d]" "${DEP_CS}" "$(($(echo "$splitArray" | grep -o ',' | wc -l) + 1 ))"
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
	mergeDepFile=$(find . -type f -iwholename "*/mergeDeps/merge_*.sh")
	for file in $mergeDepFiles; do
		mergeDeps=$(appendList "$mergeDeps" "$(cat $file)" ":")
		rm $file
	done
	numDeps=$(echo $mergeDeps | wc -w)
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
		
		printf "%d:%-10s " $i "$contig"
		
		# Gather merge inputs
		mergeOutput=merged/${contig}/merged.bam
		 markOutput=markdup/${contig}/markdupped.bam
		 baseOutput=baserecal/${contig}/baserecal.firstpass
		printOutput=printreads/${contig}/printreads.bam
		  catInputs=$(appendList "$catInputs" "${printOutput}" " ")
		depthOutput=depth/${contig}/depth #.sample_summary
		haploOutput=haplo/${contig}/${IDN}.g.vcf.gz
		
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
	done
	
	mergeReadCount=$(echo ${catInputs} | wc -w)
	
	printf "\nSM: Merge [%s] " $mergeArray
	
	if [ "$mergeArray" != "" ]; then
		DEP_CM=$(sbatch -J MC_${IDN} -a $mergeArray $(depCheck ${mergeDeps}) ${SLSBIN}/mergecontigs.sl | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CM" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_CM}
		fi
	else
		printf "done "
	fi
	
	printf "\n"
	
	exit 0
	
	printf "%s " "-> MarkDup"
	
	if [ "$markArray" != "" ]; then
		DEP_MD=$(sbatch -J MD_${IDN} -a $markArray $(depCheck ${DEP_CM}) ${SLSBIN}/markduplicates.sl | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_MD" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_MD}
		fi
	else
		printf "done "
	fi
	
	printf "%s " "-> BaseRecal"
	
	if [ "$baseArray" != "" ]; then
		DEP_BR=$(sbatch -J BR_${IDN} -a $baseArray $(depCheck ${DEP_MD}) ${SLSBIN}/baserecalibrator.sl ${markDupliOutput} ${contig} ${baseRecalOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_BR" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_BR}
		fi
	else
		printf "done "
	fi
	
	printf "%s " "-> PrintReads"
	
	if [ "$printArray" != "" ]; then
		DEP_PR=$(sbatch -J PR_${IDN} -a $printArray $(depCheck ${DEP_BR}) ${SLSBIN}/printreads.sl ${markDupliOutput} ${contig} ${baseRecalOutput} ${printReadOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_PR" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_PR}
		fi
	else
		printf "done "
	fi
	
	printf "[%02d] " $mergeReadCount
	
	# Special cases for X and Y depth of covereage as the X/YPAR1 and X/YPAR2 regions are distorted.
	# Genomic Y is rife with repeat sequences that inflate coverage so use AV5 region for that.
	# X: █▄▄▄▄▄▄▄█
	# Y: _▄▄▄▄▄▄▄_
	if [ "${contig}" == "X" ]; then
		platformFile=${genderBED}
		actualContig=${TRUEX}
	elif [ "${contig}" == "Y" ]; then
		platformFile=${genderBED}
		actualContig=${TRUEY}
	else
		platformFile=${platformBED}
		actualContig=${contig}
	fi
	
	printf "%s " "-> Depth"
	
	if [ "$depthArray" != "" ]; then
		DEP_DC=$(sbatch -J DC_${IDN} -a $depthArray $(depCheck ${DEP_PR}) ${SLSBIN}/depthofcoverage.sl ${printReadOutput} ${actualContig} ${platformFile} ${depthVirtOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_DC" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_DC}
		fi
	else
		printf "done "
	fi
	
	printf "%s " "-> Haplo"
	
	if [ "$contig" != "X" ] && [ "$contig" != "Y" ] && [ "$contig" != "MT" ]; then
		# Isn't gender or mitochondrial contig.
		# Is autosomal or GL*
		
		if [ "$haploArray" != "" ]; then
			DEP_HC=$(sbatch -J HC_${IDN} -a $haploArray $(depCheck ${DEP_PR}) ${SLSBIN}/haplotypecaller.sl ${printReadOutput} ${contig} ${haplotypeOutput} | awk '{print $4}')
			if [ $? -ne 0 ] || [ "$DEP_HC" == "" ]; then
				printf "FAILED!\n"
				exit 1
			else
				printf "%d " ${DEP_HC}
			fi
		else
			printf "done "
		fi
	fi
	
	CatVarDeps=$(appendList "$CatVarDeps" "${DEP_HC}")
	
	printf "\n"
	
	if [ ! -e coverage.sh.done ]; then
		DEP_GENDER=$(sbatch -J CGD_${IDN} $(depCheck ${DEP_DC}) ${SLSBIN}/coverage.sl ${IDN} ${READGROUP} ${PLATFORM} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_GENDER" == "" ]; then
			echo "CB: Gender Determination FAILED!"
			exit 1
		else
			echo "CB: Gender Determination ${DEP_GENDER} deps [${genderDeps}]"
		fi
	else
		echo "CB: Gender Determination done"
	fi
	
	if [ ! -e ${IDN}.coverage.sh.transfer.done ]; then
		DEP_SGENDER=$(sbatch -J TCDG_${IDN} $(depCheck ${DEP_GENDER}) ${SLSBIN}/transfer.sl ${IDN} coverage.sh ${IDN}.coverage.sh | awk '{print $4}')
	fi
	echo "CB: Save Coverage map and Gender Determination: $(depDone ${DEP_SGENDER}). Depends on \"${DEP_GENDER}\""
	
	haploXInput=printreads/X/${IDN}.bam
	haploYInput=printreads/Y/${IDN}.bam
	
	mkdir -p haplo/X/${XPAR1}
	mkdir -p haplo/X/${TRUEX}
	mkdir -p haplo/X/${XPAR2}
	mkdir -p haplo/Y
	haploXPar1Output=haplo/X/${XPAR1}/${IDN}.g.vcf.gz
	haploTRUEXOutput=haplo/X/${TRUEX}/${IDN}.g.vcf.gz
	haploXPar2Output=haplo/X/${XPAR2}/${IDN}.g.vcf.gz
		haploYOutput=haplo/Y/${IDN}.g.vcf.gz
	
	if [ ! -e ${haploXPar1Output}.done ]; then
		DEP_HCXPAR1=$(sbatch -J HC_${IDN}_${XPAR1} $(depCheck ${DEP_GENDER}) ${SLSBIN}/haplotypecaller.sl  ${haploXInput} "${XPAR1}" ${haploXPar1Output} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HCXPAR1" == "" ]; then
			echo "CB: HaplotypeCaller for ${XPAR1} FAILED!"
			exit 1
		else
			echo "CB: HaplotypeCaller for ${XPAR1} ${DEP_HCXPAR1} depends on ${DEP_GENDER}"
		fi
	else
		echo "CB: HaplotypeCaller for ${XPAR1} done"
	fi
	
	if [ ! -e ${haploTRUEXOutput}.done ]; then
		DEP_HCTRUEX=$(sbatch -J HC_${IDN}_${TRUEX} $(depCheck ${DEP_GENDER}) ${SLSBIN}/haplotypecaller.sl  ${haploXInput} "${TRUEX}" ${haploTRUEXOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HCTRUEX" == "" ]; then
			echo "CB: HaplotypeCaller for ${TRUEX} FAILED!"
			exit 1
		else
			echo "CB: HaplotypeCaller for ${TRUEX} ${DEP_HCTRUEX} depends on ${DEP_GENDER}"
		fi
	else
		echo "CB: HaplotypeCaller for ${TRUEX} done"
	fi
	
	if [ ! -e ${haploXPar2Output}.done ]; then
		DEP_HCXPAR2=$(sbatch -J HC_${IDN}_${XPAR2} $(depCheck ${DEP_GENDER}) ${SLSBIN}/haplotypecaller.sl  ${haploXInput} "${XPAR2}" ${haploXPar2Output} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HCXPAR2" == "" ]; then
			echo "CB: HaplotypeCaller for ${XPAR2} FAILED!"
			exit 1
		else
			echo "CB: HaplotypeCaller for ${XPAR2} ${DEP_HCXPAR2} depends on ${DEP_GENDER}"
		fi
	else
		echo "CB: HaplotypeCaller for ${XPAR2} done"
	fi
	
	if [ ! -e ${haploYOutput}.done ]; then
		DEP_HCY=$(sbatch -J HC_${IDN}_Y $(depCheck ${DEP_GENDER}) ${SLSBIN}/haplotypecaller.sl  ${haploYInput} "Y" ${haploYOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_HCY" == "" ]; then
			echo "CB: HaplotypeCaller for Y FAILED!"
			exit 1
		else
			echo "CB: HaplotypeCaller for Y ${DEP_HCY} depends on ${DEP_GENDER}"
		fi
	else
		echo "CB: HaplotypeCaller for Y done"
	fi
	
	# If something went wrong we want to pick up where we left off without stray : characters.
	CatVarDeps="${DEP_HC} ${DEP_HCXPAR1} ${DEP_HCTRUEX} ${DEP_HCXPAR2} ${DEP_HCY}"
	CarVarDeps2=""
	for dep in ${CatVarDeps}; do
		if [ "$CarVarDeps2" == "" ]; then
			# Initial entry.
			CarVarDeps2=${dep}
		else
			CarVarDeps2="${CarVarDeps2}:${dep}"
		fi
	done
	CatVarDeps="${CarVarDeps2}"
	
	CatVarInputs=""
	for contig in ${CONTIGS}; do
		if [ "$contig" == "MT" ]; then
			continue
		elif [ "$contig" == "X" ]; then
			thisInput="haplo/X/${XPAR1}/${IDN}.g.vcf.gz haplo/X/${TRUEX}/${IDN}.g.vcf.gz haplo/X/${XPAR2}/${IDN}.g.vcf.gz"
		else
			thisInput="haplo/${contig}/${IDN}.g.vcf.gz"
		fi
		
		if [ "$CatVarInputs" == "" ]; then
			CatVarInputs="${thisInput}"
		else
			CatVarInputs="${CatVarInputs} ${thisInput}"
		fi
	done
	
	mergePrintOutput=${IDN}.bam
	catvarOutput=${IDN}.g.vcf.gz
	
	# Merge print-read bams.
	if [ ! -e ${mergePrintOutput}.done ]; then
		DEP_CPR=$(sbatch -J CPR_${IDN} $(depCheck ${DEP_PR}) ${SLSBIN}/catreads.sl "${catInputs}" ${mergePrintOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CPR" == "" ]; then
			echo "CB: Cat Reads FAILED!"
			exit 1
		else
			echo "CB: Cat Reads  ${DEP_CPR} depends on ${mergeReadDeps} [$(echo $mergeReadDeps | sed -e 's/:/ /g' | wc -w)]"
		fi
	else
		echo "CB: Cat Reads  done"
	fi
	
	if [ ! -e ${mergePrintOutput%.bam}.bai.done ]; then
		DEP_CPRI=$(sbatch -J CPRI_${IDN} $(depCheck ${DEP_CPR}) ${SLSBIN}/catreadsindex.sl ${mergePrintOutput} ${mergePrintOutput%.bam}.bai | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CPRI" == "" ]; then
			echo "CB: Cat Reads Index FAILED!"
			exit 1
		else
			echo "CB: Cat Reads Index ${DEP_CPRI} depends on ${DEP_CPR}"
		fi
	else
		echo "CB: Cat Reads Index done"
	fi
	
	if [ ! -e ${mergePrintOutput}.transfer.done ]; then
		DEP_TCPR=$(sbatch -J TCPR_${IDN} $(depCheck ${DEP_CPR}) ${SLSBIN}/transfer.sl ${IDN} ${mergePrintOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_TCPR" == "" ]; then
			echo "CB: Save Cat Reads FAILED!"
			exit 1
		else
			echo "CB: Save Cat Reads ${DEP_TCPR} depends on ${DEP_CPR}"
		fi
	else
		echo "CB: Save Cat Reads done"
	fi
	
	if [ ! -e ${mergePrintOutput%.bam}.bai.transfer.done ]; then
		DEP_TMRI=$(sbatch -J TCPRI_${IDN} $(depCheck ${DEP_CPRI}) ${SLSBIN}/transfer.sl ${IDN} ${mergePrintOutput%.bam}.bai | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_TMRI" == "" ]; then
			echo "CB: Save Cat Reads Index FAILED!"
			exit 1
		else
			echo "CB: Save Cat Reads Index ${DEP_TMRI} depends on ${DEP_CPRI}"
		fi
	else
		echo "CB: Save Cat Reads Index done"
	fi
	
	printf "CB: Cat Variants "
	
	if [ ! -e ${catvarOutput}.done ]; then
		DEP_CV=$(sbatch -J CV_${IDN} $(depCheck ${CatVarDeps}) ${SLSBIN}/catvar.sl "${CatVarInputs}" ${catvarOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_CV" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d\n" ${DEP_CV}
		fi
	else
		printf "done\n"
	fi
	
	printf "%s " "-> Save Cat Variants " 
	
	if [ ! -e ${catvarOutput}.transfer.done ]; then
		DEP_TVF=$(sbatch -J TCV_${IDN} $(depCheck ${DEP_CV}) ${SLSBIN}/transfer.sl ${IDN} ${catvarOutput} | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_TVF" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_TVF}
		fi
	else
		printf "done "
	fi
	
	printf "& Save Variants Index "
	
	if [ ! -e ${catvarOutput}.tbi.transfer.done ]; then
		DEP_TVI=$(sbatch -J TCVI_${IDN} $(depCheck ${DEP_CV}) ${SLSBIN}/transfer.sl ${IDN} ${catvarOutput}.tbi | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_TVI" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d " ${DEP_TVI}
		fi
	else
		printf "done\n"
	fi
	
	if [ "$DEP_CV" != "" ] && [ "$DEP_CPRI" != "" ]; then
		saveDeps="${DEP_CV}:${DEP_CPRI}"
	elif [ "$DEP_CPR" != "" ]; then
		saveDeps="${DEP_CPRI}"
	elif [ "$DEP_CV" != "" ]; then
		saveDeps="${DEP_CV}"
	fi
	
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
	
	printf "%s " "-> Save commands "
	
	if [ ! -e ${IDN}.commands.txt.transfer.done ]; then
		DEP_TCMD=$(sbatch -J TCMD_${IDN} $(depCheck ${saveDeps}) ${SLSBIN}/transfer.sl ${IDN} commands.txt ${IDN}.commands.txt | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_TCMD" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%d\n" ${DEP_TCMD}
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
