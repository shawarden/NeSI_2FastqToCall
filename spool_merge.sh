#!/bin/bash

#####################################
# Align, sort and split individual  #
#  blocks by contig.                #
# Merge contigs and mark duplicates #
# Base Quality Score Recalibration  #
#####################################


   SAMPLE=${1}
 PLATFORM=${2}
MULTI_RUN=${3}

IDN=$(echo ${SAMPLE} | awk -F'[[:blank:]_]' '{print $1}')

# Get base values
source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

alignArray=""
readBlocks=$(find ./blocks -type f -iname "R1_*.fastq.gz" | wc -l)
#	purgeList="$(($readBlocks+1))-$FASTQ_MAXJOBZ"
#	scancel ${DEP_BA}_[${purgeList}] && echo "Purged excess alignment and sort jobs $purgeList"
for i in $(seq 1 ${readBlocks}); do
	if [ ! -e split/$(printf "%0${FASTQ_MAXZPAD}d" $i)/contig_split.done ]; then
		# This contig block hasn't been split yet.
		alignArray=$(appendList "$alignArray" $i ",")
	fi
done

# Dispatch alignemnt array.
# Alignemnt array doesn't have a dependency yet since ReadSplit needs to know what the aligner's JobID is to update it. Delay the start.
printf "%-22s" "Align->Sort->Split"

if [ "$alignArray" != "" ]; then
	DEP_BA=$(sbatch $(dispatch "BA") -J BA_${SAMPLE} --array=$alignArray $SLSBIN/blockalign.sl $SAMPLE | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_BA" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%sx%-2d [%s]\n" "${DEP_BA}" $(splitByChar "$alignArray" "," | wc -w) $(condenseList "$alignArray")
		echo $DEP_BA > ../lastJob.txt
	fi
else
	printf "done\n"
fi

if [ "${MULTI_RUN}" != "" ]
then
	echo "Mutliple run mode. Stopping at alignment."
	echo "Launch last run without multi-run mode enabled to complete."
	exit 0
fi

cd ..

######################################
#INJECT MULTIPLE RUNS PER INDIVIDUAL #
######################################

# This only applied to Chr X and Y and their PAR1/2 regions.
genderDeps=""

# Things CatVariants is dependant on.
# All HaplotypeCaller jobs.
CatVarDeps=""
CatVarInputs=""

# List of incomplete jobs.
mergeMarkArray=""
recalArray=""
depthArray=""
haploArray=""

# Loop though number of contigs in reference sequence.
# Build list of incomplete merged contigs.
for i in $(seq 1 ${NUMCONTIG_BLOCKS}); do
	# Build input/output file names
	contig=${CONTIGBLOCKS[$i]}	# Does bash do array lookups every time too?
	printf "%04d %-22s " $i $contig
	
	mergeMarkOutput=markdup/${contig}.bam
	    recalOutput=printreads/${contig}.bam
	  catReadInputs=$(appendList "$catReadInputs" "${recalOutput}" " ")
	    depthOutput=depth/${contig} #.sample_summary
	    haploOutput=haplo/${contig}.g.vcf.gz
	
	mkdir -p $(dirname $mergeMarkOutput)
	mkdir -p $(dirname $recalOutput)
	mkdir -p $(dirname $depthOutput)
	mkdir -p $(dirname $haploOutput)
	
	if [ ! -e ${mergeMarkOutput}.done ]; then
		mergeMarkArray=$(appendList "$mergeMarkArray"  $i ",")
		printf "MD "
	fi
	
	if [ ! -e ${recalOutput}.done ]; then
		recalArray=$(appendList "$recalArray" $i ",")
		printf "PR "
	fi
	
	if [ "$contig" != "MT" ] && [ "$contig" != "hs37d5" ] && [ "$contig" != "NC_007605" ]; then	# skip Mito, and Decoy contigs
		if [[ $contig != GL* ]]; then	# skip GL contigs for depth.
			if [ ! -e ${depthOutput}.done ]; then
				depthArray=$(appendList "$depthArray" $i ",")
				printf "DC "
			fi
		fi
		if [ "$contig" != "X" ] && [ "$contig" != "Y" ]; then	#Skip gender for haplotype.
			if [ ! -e ${haploOutput}.done ]; then
				haploArray=$(appendList "$haploArray" $i ",")
				printf "HC "
			fi
		fi
	fi
	printf "\n"
done

mergeReadCount=$(echo ${catReadInputs} | wc -w)

printf "%-22s" "Merge and Mark" 

if [ "$mergeMarkArray" != "" ]; then
	DEP_MM=$(sbatch $(dispatch "MM") -J MM_${IDN} --array $mergeMarkArray $(depCheck $DEP_BA) $SLSBIN/mergeandmark.sl | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_MM" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%sx%-2d [%s]\n" "$DEP_MM" $(splitByChar "$mergeMarkArray" "," | wc -w) $(condenseList "$mergeMarkArray")
	fi
else
	printf "done\n"
fi

printf "%-22s" "Recalibration"

if [ "$recalArray" != "" ]; then
	DEP_RC=$(sbatch $(dispatch "RC") -J RC_${IDN} --array $recalArray $(depCheckArrack $DEP_MM) $SLSBIN/recalibration.sl | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_RC" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$recalArray" "$DEP_RC" "$mergeMarkArray" "$DEP_MM"
		printf "%sx%-2d [%s]\n" "$DEP_RC" $(splitByChar "$recalArray" "," | wc -w) $(condenseList "$recalArray")
	fi
else
	printf "done\n"
fi

printf "%-22s" "Depth of Coverage"

if [ "$depthArray" != "" ]; then
	DEP_DC=$(sbatch $(dispatch "DC") -J DC_${IDN} --array $depthArray $(depCheckArrack $DEP_RC) $SLSBIN/depthofcoverage.sl $PLATFORM | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_DC" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$depthArray" "$DEP_DC" "$recalArray" "$DEP_RC"
		printf "%sx%-2d [%s]\n" "$DEP_DC" $(splitByChar "$depthArray" "," | wc -w) $(condenseList "$depthArray")
	fi
else
	printf "done\n"
fi

printf "%-22s" "HaplotypeCaller"

if [ "$haploArray" != "" ]; then
	DEP_HC=$(sbatch $(dispatch "HC") -J HC_${IDN} --array $haploArray $(depCheckArrack $DEP_RC) $SLSBIN/haplotypecaller.sl | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_HC" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		# Tie each task to the matching task in the previous array.
		tieTaskDeps "$haploArray" "$DEP_HC" "$recalArray" "$DEP_RC"
		printf "%sx%-2d [%s]\n" "$DEP_HC" $(splitByChar "$haploArray" "," | wc -w) $(condenseList "$haploArray")
		
	fi
else
	printf "done\n"
fi

printf "%-22s" "Gender Determination"

if [ ! -e coverage.sh.done ]; then
	DEP_GD=$(sbatch $(dispatch "GD") -J GD_${IDN} $(depCheck $DEP_DC) $SLSBIN/coverage.sl $IDN $PLATFORM | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_GD" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "$DEP_GD"
	fi
else
	printf "done\n"
fi

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
	DEP_HCXPAR1=$(sbatch $(dispatch "HC") -J HC_${IDN}_XPAR1 --array=23 $(depCheck $DEP_GD) $SLSBIN/haplotypecaller.sl  "$XPAR1" | awk '{print $4}')
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
	DEP_HCTRUEX=$(sbatch $(dispatch "HC") -J HC_${IDN}_TRUEX --array=23 $(depCheck $DEP_GD) $SLSBIN/haplotypecaller.sl "$TRUEX" | awk '{print $4}')
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
	DEP_HCXPAR2=$(sbatch $(dispatch "HC") -J HC_${IDN}_XPAR2 --array=23 $(depCheck $DEP_GD) $SLSBIN/haplotypecaller.sl "$XPAR2" | awk '{print $4}')
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
	DEP_HCY=$(sbatch $(dispatch "HC") -J HC_${IDN}_Y --array=24 $(depCheck $DEP_GD) ${SLSBIN}/haplotypecaller.sl "Y" | awk '{print $4}')
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
for contig in ${CONTIGBLOCKS[@]}
do
	if [ "$contig" == "MT" ] || [ "$contig" == "hs37d5" ] || [ "$contig" == "NC_007605" ]
	then	# Skip mito and decoy regions
		continue
	elif [ "$contig" == "X" ]
	then
		thisInput="haplo/${XPAR1}.g.vcf.gz haplo/${TRUEX}.g.vcf.gz haplo/${XPAR2}.g.vcf.gz"
	else
		thisInput="haplo/${contig}.g.vcf.gz"
	fi
	
	CatVarInputs=$(appendList "$CatVarInputs" "$thisInput")
done

catReadsOutput=${IDN}.bam
catVarOutput=${IDN}.g.vcf.gz

printf "%-22s" "CatReads"

# Merge print-read bams.
if [ ! -e ${catReadsOutput}.done ]; then
	DEP_CR=$(sbatch $(dispatch "CR") -J CR_${IDN} $(depCheck $DEP_RC) $SLSBIN/catreads.sl "$catReadInputs" $catReadsOutput | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_CR" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "$DEP_CR"
	fi
else
	printf "done\n"
	
	printf "%-22s" "Reads Index"
	
	if [ ! -e ${catReadsOutput%.bam}.bai.done ]; then
		DEP_RI=$(sbatch $(dispatch "RI") -J RI_${IDN} $(depCheck $DEP_CR) $SLSBIN/catreadsindex.sl $catReadsOutput ${catReadsOutput%.bam}.bai | awk '{print $4}')
		if [ $? -ne 0 ] || [ "$DEP_RI" == "" ]; then
			printf "FAILED!\n"
			exit 1
		else
			printf "%s\n" "$DEP_RI"
		fi
	else
		printf "done\n"
	fi
fi

# Add cleanup dependency.
saveDeps=$(appendList "$saveDeps" "$DEP_CR" ":")

printf "%-22s" "CatVariants"

if [ ! -e ${catVarOutput}.done ]; then
	DEP_CV=$(sbatch $(dispatch "CV") -J CV_${IDN} $(depCheck $CatVarDeps) $SLSBIN/catvar.sl "$CatVarInputs" $catVarOutput | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_CV" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "$DEP_CV"
	fi
else
	printf "done\n" 
fi

# Add cleanup dependency.
saveDeps=$(appendList "$saveDeps" "$DEP_CV" ":")

printf "%-22s" "Save workflow"

if [ "$saveDeps" != "" ]; then
	DEP_CU=$(sbatch $(dispatch "TF") -J CU_${IDN} $(depCheck $saveDeps) $SLSBIN/clearup.sl $IDN | awk '{print $4}')
	if [ $? -ne 0 ] || [ "$DEP_CU" == "" ]; then
		printf "FAILED!\n"
		exit 1
	else
		printf "%s\n" "$DEP_CU"
		echo $DEP_CU > lastJob.txt
	fi
else
	printf "done\n"
fi
