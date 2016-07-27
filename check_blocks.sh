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

alignOutput=aligned/${READGROUP}_${BLOCK}.bam
sortOutput=sorted/${READGROUP}_${BLOCK}.bam

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

function batchJob {
	if [ ! -e ${alignOutput}.done ]; then
		DEP_PA=$(sbatch -J PA_${SAMPLE}_${BLOCK} ${SLSBIN}/align.sl ${SAMPLE} ${READGROUP} ${BLOCK} ${alignOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Alignment FAILED!"
			exit 1
		else
			echo "CheckBlock: Alignment ${DEP_PA}"
		fi
	else
		echo "CheckBlock: Alignment already COMLETED!"
	fi
	
	if [ ! -e ${sortOutput}.done ]; then
		DEP_SS=$(sbatch -J SS_${SAMPLE}_${BLOCK} $(depCheck ${DEP_PA}) ${SLSBIN}/sortsam.sl ${SAMPLE} ${alignOutput} ${sortOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: SortSam FAILED!"
			exit 1
		else
			echo "CheckBlock: SortSam ${DEP_SS} depends on ${DEP_PA}"
		fi
	else
		echo "CheckBlock: SortSam already COMLETED!"
	fi
	
	printf "CheckBlock: Contig split "
	
	for contig in ${CONTIGS}; do
		# Load existing dependency list for this contig.
		if [ -e mergeDeps/${contig}.sh ]; then
			source mergeDeps/${contig}.sh
		fi
		
		mkdir -p contig/${contig}
		
		splitOutput=contig/${contig}/${READGROUP}_${BLOCK}.bam
		
		if [ ! -e ${splitOutput}.done ]; then
			DEP_CS=$(sbatch -J CS_${SAMPLE}_${BLOCK}_${contig} $(depCheck ${DEP_SS}) ${SLSBIN}/contigsplit.sl ${SAMPLE} ${contig} ${sortOutput} ${splitOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				printf "%10s FAILED!\n" "${contig} "
				exit 1
			else
				printf "%10s %s depends on %s" "${contig}" "${DEP_CS}" "${DEP_SS} "
			fi
			
			# Add to list of dependencies for this contig
			if [ "$DEP_CS" != "" ]; then
				if [ "$mergeDeps" == "" ]; then
					# Initial entry.
					mergeDeps="${DEP_CS}"
				else
					# Additional entry.
					mergeDeps="${mergeDeps}:${DEP_CS}"
				fi
			fi
			
			# Store list of dependencies for this contig.
			mkdir -p mergeDeps
			echo "mergeDeps=${mergeDeps}" > mergeDeps/${contig}.sh
		else
			printf "%10s already COMLETED!" "${contig}"
		fi
	done
	echo ""
	
	if [ "$BLOCK" == "$NEXT" ]; then	# This is the final chunk so we can spool up the next step
		echo "CheckBlock: Last split so spool up block merging."
		spoolmerge
	fi
}

function spoolmerge {
	READGROUP=$(cat blocks/R1_ReadGroup.txt)
	
	cd ../
	
	mkdir -p slurm
	
	if [ $(echo ${READGROUP} | wc -w) -gt 1 ]
	then
		echo "CheckBlock: Too many read-groups. No spooling further pipeline."
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
	
	for contig in ${CONTIGS}; do
		mergeDepFile=$(find . -type f -iwholename "*/mergeDeps/${contig}.sh")
		echo "CheckBlock: mergeDepFiles for ${contig}: [${mergeDepFile}]"
		if [ "$mergeDepFile" != "" ]; then
			source ${mergeDepFile}
			rm ${mergeDepFile}
		fi
		
		# Reset dependencies
		DEP_CM=""
		DEP_MD=""
		DEP_BR=""
		DEP_PR=""
		DEP_DC=""
		DEP_HC=""
		
		numDeps=$(echo $mergeDeps | wc -w)
		numContigs=$(echo $CONTIGS | wc -w)
		echo "CheckBlock: Dependencies for ${contig}: $numDeps Contigs: $numContigs"
		
		mkdir -p merged/${contig}
		mergeContOutput=merged/${contig}/${IDN}.bam
		
		mkdir -p markdup/${contig}
		markDupliOutput=markdup/${contig}/${IDN}.bam
		
		mkdir -p baserecal/${contig}
		baseRecalOutput=baserecal/${contig}/${IDN}.firstpass
		
		mkdir -p printreads/${contig}
		printReadOutput=printreads/${contig}/${IDN}.bam
		
		mkdir -p depth/${contig}
		depthVirtOutput=depth/${contig}/${IDN} #.sample_summary
		depthActuOutput=${depthOutput}.sample_summary
		
		mkdir -p haplo/${contig}
		haplotypeOutput=haplo/${contig}/${IDN}.g.vcf.gz
		
		if [ ! -e ${mergeContOutput}.done ]; then
			DEP_CM=$(sbatch -J MC_${IDN}_${contig} $(depCheck ${mergeDeps}) ${SLSBIN}/mergecontigs.sl ${contig} ${mergeContOutput} "yes" | awk '{print $4}')
			if [ $? -ne 0 ]; then
				echo "CheckBlock: Merge Contig for ${contig} FAILED!"
				exit 1
			else
				echo "CheckBlock: Merge Contig for ${contig} ${DEP_CM} depends on ${mergeDeps}"
			fi
		else
			echo "CheckBlock: Merge Contig for ${contig} already COMLETED!"
		fi
		
		if [ ! -e ${markDupliOutput}.done ]; then
			DEP_MD=$(sbatch -J MD_${IDN}_${contig} $(depCheck ${DEP_CM}) ${SLSBIN}/markduplicates.sl ${mergeContOutput} ${markDupliOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				echo "CheckBlock: Mark Duplicates for ${contig} FAILED!"
				exit 1
			else
				echo "CheckBlock: Mark Duplicates for ${contig} ${DEP_MD} depends on ${DEP_CM}"
			fi
		else
			echo "CheckBlock: Mark Duplicates for ${contig} already COMLETED!"
		fi
		
		if [ ! -e ${baseRecalOutput}.done ]; then
			DEP_BR=$(sbatch -J BR_${IDN}_${contig} $(depCheck ${DEP_MD}) ${SLSBIN}/baserecalibrator.sl ${markDupliOutput} ${contig} ${baseRecalOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				echo "CheckBlock: Base Recalibration for ${contig} FAILED!"
				exit 1
			else
				echo "CheckBlock: Base Recalibration for ${contig} ${DEP_BR} depends on ${DEP_MD}"
			fi
		else
			echo "CheckBlock: Base Recalibration for ${contig} already COMLETED!"
		fi
		
		if [ ! -e ${printReadOutput}.done ]; then
			DEP_PR=$(sbatch -J PR_${IDN}_${contig} $(depCheck ${DEP_BR}) ${SLSBIN}/printreads.sl ${markDupliOutput} ${contig} ${baseRecalOutput} ${printReadOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				echo "CheckBlock: Print Reads for ${contig} FAILED!"
				exit 1
			else
				echo "CheckBlock: Print Reads for ${contig} ${DEP_PR} depends on ${DEP_BR}"
			fi
		else
			echo "CheckBlock: Print Reads for ${contig} already COMLETED!"
		fi
		
		# Build list of print read merge inputs and dependencies
		if [ "$mergeReadInputs" == "" ]; then
			# Initial entry.
			mergeReadInputs="${printReadOutput}"
		else
			# Additional entry.
			mergeReadInputs="${mergeReadInputs} ${printReadOutput}"
		fi
		# No quotes so no newline chars
		echo "CheckBlock: Merge inputs for ${contig}: $(echo ${mergeReadInputs} | wc -w)"
		
		if [ "$DEP_PR" != "" ]
		then
			if [ "$mergeReadDeps" == "" ]
			then
				[ "$DEP_PR" != "" ] && mergeReadDeps="$DEP_PR"
			else
				[ "$DEP_PR" != "" ] && mergeReadDeps="$mergeReadDeps:$DEP_PR"
			fi
			echo "CheckBlock: Merge dependencies for ${contig}: \"${mergeReadDeps}\""
		else
			echo "CheckBlock: No merge dependencies for ${contig} as Print-reads already run."
		fi
		
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
		
		if [ ! -e ${depthVirtOutput}.done ]; then
			DEP_DC=$(sbatch -J DC_${IDN}_${actualContig} $(depCheck ${DEP_PR}) ${SLSBIN}/depthofcoverage.sl ${printReadOutput} ${actualContig} ${platformFile} ${depthVirtOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				echo "CheckBlock: Depth of Coverage for ${contig} FAILED!"
				exit 1
			else
				echo "CheckBlock: Depth of Coverage for ${contig} ${DEP_DC} depends on ${DEP_PR}"
			fi
			
			if [ "$contig" != "MT" ] && [ "$DEP_DC" != "" ]; then
				# This is a non-gl and non-mt contig.
				# Gender determination depends on this.
				
				# Build list of gender determination dependencies.
				if [ "$genderDeps" == "" ]; then
					# Initial entry.
					genderDeps="${DEP_DC}"
				else
					# Additional entry.
					genderDeps="${genderDeps}:${DEP_DC}"
				fi
			fi
		else
			echo "CheckBlock: Depth of Coverage for ${contig} already COMLETED!"
		fi
		
		if [ "$contig" != "X" ] && [ "$contig" != "Y" ] && [ "$contig" != "MT" ]; then
			# Isn't gender or mitochondrial contig.
			# Is autosomal or GL*
			
			if [ ! -e ${haplotypeOutput}.done ]; then
				DEP_HC=$(sbatch -J HC_${IDN}_${contig} $(depCheck ${DEP_PR}) ${SLSBIN}/haplotypecaller.sl ${printReadOutput} ${contig} ${haplotypeOutput} | awk '{print $4}')
				if [ $? -ne 0 ]; then
					echo "CheckBlock: HaplotypeCaller for ${contig} FAILED!"
					exit 1
				else
					echo "CheckBlock: HaplotypeCaller for ${contig} ${DEP_HC} depends on ${DEP_PR}"
				fi
			else
				echo "CheckBlock: HaplotypeCaller for ${contig} already COMLETED!"
			fi
		fi
		
		if [ "${DEP_HC}" != "" ]; then
			# Picked up a dependency so add it to the catvar list.
			if [ "$CatVarDeps" == "" ]; then
				# Initial entry.
				CatVarDeps="${DEP_HC}"
			else 
				# Additional entry.
				CatVarDeps="${CatVarDeps} ${DEP_HC}"
			fi
		fi
	done
	
	if [ ! -e coverage.sh.done ]; then
		DEP_GENDER=$(sbatch -J CGD_${IDN} $(depCheck ${genderDeps}) ${SLSBIN}/coverage.sl ${IDN} ${READGROUP} ${PLATFORM} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Coverage & Gender determination FAILED!"
			exit 1
		else
			echo "CheckBlock: Coverage & Gender determination ${DEP_GENDER} depends on ${genderDeps}"
		fi
	else
		echo "CheckBlock: Coverage & Gender determination already COMLETED!"
	fi
	
	if [ ! -e ${IDN}.coverage.sh.transfer.done ]; then
		DEP_SGENDER=$(sbatch -J TCDG_${IDN} $(depCheck ${DEP_GENDER}) ${SLSBIN}/transfer.sl ${IDN} coverage.sh ${IDN}.coverage.sh | awk '{print $4}')
	fi
	echo "CheckBlock: Save Coverage map and Gender Determination: $(depDone ${DEP_SGENDER}). Depends on \"${DEP_GENDER}\""
	
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
		if [ $? -ne 0 ]; then
			echo "CheckBlock: HaplotypeCaller for ${XPAR1} FAILED!"
			exit 1
		else
			echo "CheckBlock: HaplotypeCaller for ${XPAR1} ${DEP_HCXPAR1} depends on ${DEP_GENDER}"
		fi
	else
		echo "CheckBlock: HaplotypeCaller for ${XPAR1} already COMLETED!"
	fi
	
	if [ ! -e ${haploTRUEXOutput}.done ]; then
		DEP_HCTRUEX=$(sbatch -J HC_${IDN}_${TRUEX} $(depCheck ${DEP_GENDER}) ${SLSBIN}/haplotypecaller.sl  ${haploXInput} "${TRUEX}" ${haploTRUEXOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: HaplotypeCaller for ${TRUEX} FAILED!"
			exit 1
		else
			echo "CheckBlock: HaplotypeCaller for ${TRUEX} ${DEP_HCTRUEX} depends on ${DEP_GENDER}"
		fi
	else
		echo "CheckBlock: HaplotypeCaller for ${TRUEX} already COMLETED!"
	fi
	
	if [ ! -e ${haploXPar2Output}.done ]; then
		DEP_HCXPAR2=$(sbatch -J HC_${IDN}_${XPAR2} $(depCheck ${DEP_GENDER}) ${SLSBIN}/haplotypecaller.sl  ${haploXInput} "${XPAR2}" ${haploXPar2Output} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: HaplotypeCaller for ${XPAR2} FAILED!"
			exit 1
		else
			echo "CheckBlock: HaplotypeCaller for ${XPAR2} ${DEP_HCXPAR2} depends on ${DEP_GENDER}"
		fi
	else
		echo "CheckBlock: HaplotypeCaller for ${XPAR2} already COMLETED!"
	fi
	
	if [ ! -e ${haploYOutput}.done ]; then
		DEP_HCY=$(sbatch -J HC_${IDN}_Y $(depCheck ${DEP_GENDER}) ${SLSBIN}/haplotypecaller.sl  ${haploYInput} "Y" ${haploYOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: HaplotypeCaller for Y FAILED!"
			exit 1
		else
			echo "CheckBlock: HaplotypeCaller for Y ${DEP_HCY} depends on ${DEP_GENDER}"
		fi
	else
		echo "CheckBlock: HaplotypeCaller for Y already COMLETED!"
	fi
	
	# If something went wrong we want to pick up where we left off without stray : characters.
	CatVarDeps="${CatVarDeps} ${DEP_HCXPAR1} ${DEP_HCTRUEX} ${DEP_HCXPAR2} ${DEP_HCY}"
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
#	if [ ! -e ${mergePrintOutput}.done ]; then
#		DEP_MR=$(sbatch -J MR_${IDN} $(depCheck ${mergeReadDeps}) ${SLSBIN}/mergereads.sl "${mergeReadInputs}" ${mergePrintOutput} | awk '{print $4}')
#	fi
#	echo "CheckBlock: Merge Reads $(depDone ${DEP_MR}). Depends on \"${mergeReadDeps}\""
	
	if [ ! -e ${mergePrintOutput}.done ]; then
		DEP_CPR=$(sbatch -J CPR_${IDN} $(depCheck ${mergeReadDeps}) ${SLSBIN}/catreads.sl "${mergeReadInputs}" ${mergePrintOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Cat Reads FAILED!"
			exit 1
		else
			echo "CheckBlock: Cat Reads  ${DEP_CPR} depends on ${mergeReadDeps}"
		fi
	else
		echo "CheckBlock: Cat Reads  already COMLETED!"
	fi
	
	if [ ! -e ${mergePrintOutput%.bam}.bai.done ]; then
		DEP_CPRI=$(sbatch -J CPRI_${IDN} $(depCheck ${DEP_CPR}) ${SLSBIN}/catreadsindex.sl ${mergePrintOutput} ${mergePrintOutput%.bam}.bai | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Cat Reads Index FAILED!"
			exit 1
		else
			echo "CheckBlock: Cat Reads Index ${DEP_CPRI} depends on ${DEP_CPR}"
		fi
	else
		echo "CheckBlock: Cat Reads Index already COMLETED!"
	fi
	
	if [ ! -e ${mergePrintOutput}.transfer.done ]; then
		DEP_TCPR=$(sbatch -J TCPR_${IDN} $(depCheck ${DEP_CPR}) ${SLSBIN}/transfer.sl ${IDN} ${mergePrintOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Save Cat Reads FAILED!"
			exit 1
		else
			echo "CheckBlock: Save Cat Reads ${DEP_TCPR} depends on ${DEP_CPR}"
		fi
	else
		echo "CheckBlock: Save Cat Reads already COMLETED!"
	fi
	
	if [ ! -e ${mergePrintOutput%.bam}.bai.transfer.done ]; then
		DEP_TMRI=$(sbatch -J TCPRI_${IDN} $(depCheck ${DEP_CPRI}) ${SLSBIN}/transfer.sl ${IDN} ${mergePrintOutput%.bam}.bai | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Save Cat Reads Index FAILED!"
			exit 1
		else
			echo "CheckBlock: Save Cat Reads Index ${DEP_TMRI} depends on ${DEP_CPRI}"
		fi
	else
		echo "CheckBlock: Save Cat Reads Index already COMLETED!"
	fi
	
	if [ ! -e ${catvarOutput}.done ]; then
		DEP_CV=$(sbatch -J CV_${IDN} $(depCheck ${CatVarDeps}) ${SLSBIN}/catvar.sl "${CatVarInputs}" ${catvarOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Cat Variants FAILED!"
			exit 1
		else
			echo "CheckBlock: Cat Variants ${DEP_CV} depends on ${CatVarDeps}"
		fi
	else
		echo "CheckBlock: Cat Variants already COMLETED!"
	fi
	
	if [ ! -e ${catvarOutput}.transfer.done ]; then
		DEP_TVF=$(sbatch -J TCV_${IDN} $(depCheck ${DEP_CV}) ${SLSBIN}/transfer.sl ${IDN} ${catvarOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Save Cat Variants FAILED!"
			exit 1
		else
			echo "CheckBlock: Save Cat Variants ${DEP_TVF} depends on ${DEP_CV}"
		fi
	else
		echo "CheckBlock: Save Cat Variants already COMLETED!"
	fi
	
	if [ ! -e ${catvarOutput}.tbi.transfer.done ]; then
		DEP_TVI=$(sbatch -J TCVI_${IDN} $(depCheck ${DEP_CV}) ${SLSBIN}/transfer.sl ${IDN} ${catvarOutput}.tbi | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Save Cat Variants Index FAILED!"
			exit 1
		else
			echo "CheckBlock: Save Cat Variants Index ${DEP_TVI} depends on ${DEP_CV}"
		fi
	else
		echo "CheckBlock: Save Cat Variants Index already COMLETED!"
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
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Save metrics FAILED!"
			exit 1
		else
			echo "CheckBlock: Save metrics ${DEP_TMet} depends on ${saveDeps}"
		fi
	else
		echo "CheckBlock: Save metrics already COMLETED!"
	fi
	
	if [ ! -e ${IDN}.commands.txt.transfer.done ]; then
		DEP_TCMD=$(sbatch -J TCMD_${IDN} $(depCheck ${saveDeps}) ${SLSBIN}/transfer.sl ${IDN} commands.txt ${IDN}.commands.txt | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CheckBlock: Save commands FAILED!"
			exit 1
		else
			echo "CheckBlock: Save commands ${DEP_TCMD} depends on ${saveDeps}"
		fi
	else
		echo "CheckBlock: Save commands already COMLETED!"
	fi
}

if [ "$READNUM" == "R1" ]; then	# BLOCK Read1 has completed!
	if [ -e ${curRead2File}.done ]; then	# NEXT Read2 exists so BLOCK Read2 has completed!
		echo "CheckBlock: R1 ${BLOCK} complete and R2 ${NEXT} exists so R2 ${BLOCK} is also complete!" #| tee -a check_${BLOCK}.txt
		batchJob
	else	# NEXT Read2 doesn't exist yet so BLOCK Read2 hasn't completed!
		echo "CheckBlock: R1 ${BLOCK} completed but R2 ${NEXT} doesn't exist yet so R2 ${BLOCK} hasn't finished!" #| tee -a check_${BLOCK}.txt
	fi
elif [ "$READNUM" == "R2" ]; then	# BLOCK Read2 has completed!
	if [ -e ${curRead1File}.done ]; then		# NEXT Read1 exists so BLOCK Read1 has completed!
		echo "CheckBlock: R2 ${BLOCK} complete and R1 ${NEXT} exists so R1 ${BLOCK} is also complete!" #| tee -a check_${BLOCK}.txt
		batchJob
	else	# NEXT Read1 doesn't exist yet so BLOCK Read1 hasn't completed!
		echo "CheckBlock: R2 ${BLOCK} complete but R1 ${NEXT} doesn't exist yet so R1 ${BLOCK} hasn't finished!" #| tee -a check_${BLOCK}.txt
	fi
fi
