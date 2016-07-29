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

function spoolAlign {
	outLine="SA: Alignment"
	
	if [ ! -e ${alignOutput}.done ]; then
		DEP_PA=$(sbatch -J PA_${SAMPLE}_${BLOCK} ${SLSBIN}/align.sl ${SAMPLE} ${READGROUP} ${BLOCK} ${alignOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			outLine=$(appendList "$outLine" "FAILED!")
			exit 1
		else
			outLine=$(appendList "$outLine" "${DEP_PA}")
		fi
	else
		outLine=$(appendList "$outLine" "done")
	fi
	
	outLine=$(appendList "$outLine" "-> SortSAM")
	
	if [ ! -e ${sortOutput}.done ]; then
		DEP_SS=$(sbatch -J SS_${SAMPLE}_${BLOCK} $(depCheck ${DEP_PA}) ${SLSBIN}/sortsam.sl ${SAMPLE} ${alignOutput} ${sortOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			outLine=$(appendList "$outLine" "FAILED!")
			exit 1
		else
			outLine=$(appendList "$outLine" "${DEP_SS}")
		fi
	else
		outLine=$(appendList "$outLine" "done")
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
	
	outLine=$(appendList "$outLine" "-> ContigSplit")
	
	if [ "$splitArray" != "" ]; then
		# Split elements defined!
		
		DEP_CS=$(sbatch -J CS_${SAMPLE}_${BLOCK} -a $splitArray $(depCheck ${DEP_SS}) ${SLSBIN}/contigsplit.sl ${SAMPLE} ${sortOutput} ${BLOCK} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			outLine=$(appendList "$outLine" "FAILED!")
			exit 1
		else
			outLine=$(appendList "$outLine" "${DEP_CS}[$(($(echo "$splitArray" | grep -o ',' | wc -l) + 1 ))]")
		fi
		mkdir -p mergeDeps
		echo "${DEP_CS}" > mergeDeps/merge_${DEP_CS}.sh
	else
		# No array elements define. They're all done!
		outLine=$(appendList "$outLine" "done!")
	fi
	
	echo "${outLine}"
	
	if [ "$BLOCK" == "$NEXT" ]; then	# This is the final chunk so we can spool up the next step
		echo "CB: Last split so spool up block merging."
		# spoolMerge
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
	echo "CB: Merge Dependencies: [${mergeDeps}]"
	
	for contig in ${CONTIGS}; do
		outLine="SM: ${contig} Merge"
		
		# Reset dependencies
		DEP_CM=""
		DEP_MD=""
		DEP_BR=""
		DEP_PR=""
		DEP_DC=""
		DEP_HC=""
		
		mergeContOutput=merged/${contig}/merged.bam
		mkdir -p $(dirname $mergeContOutput)
		
		markDupliOutput=markdup/${contig}/markdupped.bam
		mkdir -p $(dirname $markDupliOutput)
		
		baseRecalOutput=baserecal/${contig}/baserecal.firstpass
		mkdir -p $(dirname $baseRecalOutput)
		
		printReadOutput=printreads/${contig}/printreads.bam
		mkdir -p $(dirname $printReadOutput)
		
		depthVirtOutput=depth/${contig}/${IDN} #.sample_summary
		depthActuOutput=${depthVirtOutput}.sample_summary
		mkdir -p $(dirname $depthVirtOutput)
		
		haplotypeOutput=haplo/${contig}/${IDN}.g.vcf.gz
		mkdir -p $(dirname $haplotypeOutput)
		
		numDeps=$(echo $mergeDeps | wc -w)
		numContigs=$(echo $CONTIGS | wc -w)
		echo "CB: Dependencies for ${contig}: $numDeps Contigs: $numContigs"
		
		if [ ! -e ${mergeContOutput}.done ]; then
			DEP_CM=$(sbatch -J MC_${IDN}_${contig} $(depCheck ${mergeDeps}) ${SLSBIN}/mergecontigs.sl ${contig} ${mergeContOutput} "yes" | awk '{print $4}')
			if [ $? -ne 0 ]; then
				outLine=$(appendList "$outLine" "FAILED!")
				echo ${outLine}
				exit 1
			else
				outLine=$(appendList "$outLine" "${DEP_CM}")
			fi
		else
			outLine=$(appendList "$outLine" "done")
		fi
		
		outLine=$(appendList "$outLine" "MarkDup")
		
		if [ ! -e ${markDupliOutput}.done ]; then
			DEP_MD=$(sbatch -J MD_${IDN}_${contig} $(depCheck ${DEP_CM}) ${SLSBIN}/markduplicates.sl ${mergeContOutput} ${markDupliOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				outLine=$(appendList "$outLine" "FAILED!")
				echo ${outLine}
				exit 1
			else
				outLine=$(appendList "$outLine" "${DEP_MD}")
			fi
		else
			outLine=$(appendList "$outLine" "done")
		fi
		
		outLine=$(appendList "$outLine" "BaseRecal")
		
		if [ ! -e ${baseRecalOutput}.done ]; then
			DEP_BR=$(sbatch -J BR_${IDN}_${contig} $(depCheck ${DEP_MD}) ${SLSBIN}/baserecalibrator.sl ${markDupliOutput} ${contig} ${baseRecalOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				outLine=$(appendList "$outLine" "FAILED")
				echo ${outLine}
				exit 1
			else
				outLine=$(appendList "$outLine" "${DEP_BR}")
			fi
		else
			outLine=$(appendList "$outLine" "done")
		fi
		
		outLine=$(appendList "$outLine" "PrintReads")
		
		if [ ! -e ${printReadOutput}.done ]; then
			DEP_PR=$(sbatch -J PR_${IDN}_${contig} $(depCheck ${DEP_BR}) ${SLSBIN}/printreads.sl ${markDupliOutput} ${contig} ${baseRecalOutput} ${printReadOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				outLine=$(appendList "$outLine" "FAILED!")
				echo ${outLine}
				exit 1
			else
				outLine=$(appendList "$outLine" "${DEP_PR}")
			fi
		else
			outLine=$(appendList "$outLine" "done")
		fi
		
		# Build list of print read merge inputs and dependencies
		mergeReadInputs=$(appendList "$mergeReadInputs" "${printReadOutput}" " ")
		
		# No quotes so no newline chars
		outLine=$(appendList "$outLine" "merge inputs [$(echo ${mergeReadInputs} | wc -w)]")
		
		if [ "$DEP_PR" != "" ]; then
			mergeReadDeps=$(appendList "$mergeReadDeps" "$DEP_PR" ":")
			outLine=$(appendList "$outLine" "deps [${mergeReadDeps}]")
		else
			outLine=$(appendList "$outLine" "dep-less")
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
		
		outLine=$(appendList "$outLine" "Depth")
		
		if [ ! -e ${depthVirtOutput}.done ]; then
			DEP_DC=$(sbatch -J DC_${IDN}_${actualContig} $(depCheck ${DEP_PR}) ${SLSBIN}/depthofcoverage.sl ${printReadOutput} ${actualContig} ${platformFile} ${depthVirtOutput} | awk '{print $4}')
			if [ $? -ne 0 ]; then
				outLine=$(appendList "$outLine" "${contig} FAILED!")
				echo ${outLine}
				exit 1
			else
				outLine=$(appendList "$outLine" "${contig} ${DEP_DC}")
			fi
			
			if [ "$DEP_DC" != "" ]; then
				# Build list of gender determination dependencies.
				genderDeps=$(appendList "$genderDeps" "${DEP_DC}" ":")
			fi
		else
			outLine=$(appendList "$outLine" "${contig} done")
		fi
		
		outLine=$(appendList "$outLine" "Haplo")
		
		if [ "$contig" != "X" ] && [ "$contig" != "Y" ] && [ "$contig" != "MT" ]; then
			# Isn't gender or mitochondrial contig.
			# Is autosomal or GL*
			
			if [ ! -e ${haplotypeOutput}.done ]; then
				DEP_HC=$(sbatch -J HC_${IDN}_${contig} $(depCheck ${DEP_PR}) ${SLSBIN}/haplotypecaller.sl ${printReadOutput} ${contig} ${haplotypeOutput} | awk '{print $4}')
				if [ $? -ne 0 ]; then
					outLine=$(appendList "$outLine" "${contig} FAILED!")
					echo ${outLine}
					exit 1
				else
					outLine=$(appendList "$outLine" "${contig} ${DEP_HC}")
				fi
			else
				outLine=$(appendList "$outLine" "${contig} done")
			fi
		fi
		
		CatVarDeps=$(appendList "$CatVarDeps" "${DEP_HC}")
		
		echo "${outLine}"
		
	done
	
	if [ ! -e coverage.sh.done ]; then
		DEP_GENDER=$(sbatch -J CGD_${IDN} $(depCheck ${genderDeps}) ${SLSBIN}/coverage.sl ${IDN} ${READGROUP} ${PLATFORM} | awk '{print $4}')
		if [ $? -ne 0 ]; then
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
		if [ $? -ne 0 ]; then
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
		if [ $? -ne 0 ]; then
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
		if [ $? -ne 0 ]; then
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
		if [ $? -ne 0 ]; then
			echo "CB: HaplotypeCaller for Y FAILED!"
			exit 1
		else
			echo "CB: HaplotypeCaller for Y ${DEP_HCY} depends on ${DEP_GENDER}"
		fi
	else
		echo "CB: HaplotypeCaller for Y done"
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
#	echo "CB: Merge Reads $(depDone ${DEP_MR}). Depends on \"${mergeReadDeps}\""
	
	if [ ! -e ${mergePrintOutput}.done ]; then
		DEP_CPR=$(sbatch -J CPR_${IDN} $(depCheck ${mergeReadDeps}) ${SLSBIN}/catreads.sl "${mergeReadInputs}" ${mergePrintOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CB: Cat Reads FAILED!"
			exit 1
		else
			echo "CB: Cat Reads  ${DEP_CPR} depends on ${mergeReadDeps}"
		fi
	else
		echo "CB: Cat Reads  done"
	fi
	
	if [ ! -e ${mergePrintOutput%.bam}.bai.done ]; then
		DEP_CPRI=$(sbatch -J CPRI_${IDN} $(depCheck ${DEP_CPR}) ${SLSBIN}/catreadsindex.sl ${mergePrintOutput} ${mergePrintOutput%.bam}.bai | awk '{print $4}')
		if [ $? -ne 0 ]; then
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
		if [ $? -ne 0 ]; then
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
		if [ $? -ne 0 ]; then
			echo "CB: Save Cat Reads Index FAILED!"
			exit 1
		else
			echo "CB: Save Cat Reads Index ${DEP_TMRI} depends on ${DEP_CPRI}"
		fi
	else
		echo "CB: Save Cat Reads Index done"
	fi
	
	if [ ! -e ${catvarOutput}.done ]; then
		DEP_CV=$(sbatch -J CV_${IDN} $(depCheck ${CatVarDeps}) ${SLSBIN}/catvar.sl "${CatVarInputs}" ${catvarOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CB: Cat Variants FAILED!"
			exit 1
		else
			echo "CB: Cat Variants ${DEP_CV} depends on ${CatVarDeps}"
		fi
	else
		echo "CB: Cat Variants done"
	fi
	
	if [ ! -e ${catvarOutput}.transfer.done ]; then
		DEP_TVF=$(sbatch -J TCV_${IDN} $(depCheck ${DEP_CV}) ${SLSBIN}/transfer.sl ${IDN} ${catvarOutput} | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CB: Save Cat Variants FAILED!"
			exit 1
		else
			echo "CB: Save Cat Variants ${DEP_TVF} depends on ${DEP_CV}"
		fi
	else
		echo "CB: Save Cat Variants done"
	fi
	
	if [ ! -e ${catvarOutput}.tbi.transfer.done ]; then
		DEP_TVI=$(sbatch -J TCVI_${IDN} $(depCheck ${DEP_CV}) ${SLSBIN}/transfer.sl ${IDN} ${catvarOutput}.tbi | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CB: Save Cat Variants Index FAILED!"
			exit 1
		else
			echo "CB: Save Cat Variants Index ${DEP_TVI} depends on ${DEP_CV}"
		fi
	else
		echo "CB: Save Cat Variants Index done"
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
			echo "CB: Save metrics FAILED!"
			exit 1
		else
			echo "CB: Save metrics ${DEP_TMet} depends on ${saveDeps}"
		fi
	else
		echo "CB: Save metrics done"
	fi
	
	if [ ! -e ${IDN}.commands.txt.transfer.done ]; then
		DEP_TCMD=$(sbatch -J TCMD_${IDN} $(depCheck ${saveDeps}) ${SLSBIN}/transfer.sl ${IDN} commands.txt ${IDN}.commands.txt | awk '{print $4}')
		if [ $? -ne 0 ]; then
			echo "CB: Save commands FAILED!"
			exit 1
		else
			echo "CB: Save commands ${DEP_TCMD} depends on ${saveDeps}"
		fi
	else
		echo "CB: Save commands done"
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
