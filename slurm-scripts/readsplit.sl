#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		ReadSplit
#SBATCH --time			0-02:30:00
#SBATCH --mem-per-cpu	2048
#SBATCH --cpus-per-task	8
#SBATCH --constraint	avx
#SBATCH --array			1,2
#SBATCH --error			slurm/RS_%A_%a.out
#SBATCH --output		slurm/RS_%A_%a.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

 SAMPLE=${1}
READNUM=${SLURM_ARRAY_TASK_ID}
 
[ $READNUM -eq 1 ] && INPUT=${2} || INPUT=${3}

ALIGN_ARRAY=${4}
 SORT_ARRAY=${5} 

# Set entire Alignment Array dependency to this job's success.
# check_blocks script will release individual array elements, then purge the rest once BLOCK and NEXT match.
scontrol update JobID=$ALIGN_ARRAY StartTime=now Dependency=afterok:$SLURM_ARRAY_JOB_ID

HEADER="RS"

echo "$HEADER: R${READNUM} ${INPUT}"
date

# Make sure input exists!
if ! inFile; then exit $EXIT_IO; fi

module load ${MOD_ZLIB}

# Scan the FastQ file for index sequence.
# Returns the index with the highest count within FASTQ_MAXSCAN lines from top.
function getBestIndex {
	${CAT_CMD} ${INPUT} | head -${FASTQ_MAXSCAN} | awk \
	-F':' \
	'NR%4==1{words[$10]++}
	END {
		for (w in words) {
			print words[w], w
		}
	}' | sort -rn | awk '{print $2}' | head -1
}

# Splits the sourceFile into sub-files based on header and index.
# If any element of the header exect the index differs, split to new output file.
# If the index varies by more than FASTQ_MAXDIFF then split to new output file.
# Output files are appended with the index line data.
function splitByReadGroupAndCompress {
	if ! ${CAT_CMD} ${INPUT} | awk -F'[@:]' \
		-v zeroPad="$FASTQ_MAXZPAD" \
		-v outHeader="$HEADER" \
		-v sampleID="$SAMPLE" \
		-v readNumber="R$READNUM" \
		-v seqIndex="$bestIndex" \
		-v maxDelta="$FASTQ_MAXDIFF" \
		-v splitPoint="$FASTQ_MAXREAD" \
		-v compCmd="$ZIP_CMD" \
		-v alignArray="$ALIGN_ARRAY" \
		-v sortArray="$SORT_ARRAY" \
		-v pBin="$PBIN" \
		'
		BEGIN{
			curBlock=0
			blockCount=0
			padBlock=sprintf("%0"zeroPad"d", curBlock)
			if (system("[ -e blocks/*_"readNumber"_"padBlock".fastq.gz.done ]") == 0) {
				# a blocks/..._Rn_00000.fastq.gz.done file exists so skip this block
				writeBlock=0
				print outHeader": Skipping block "padBlock" as it is already completed!"
			} else {
				writeBlock=1
				print outHeader": Block "padBlock" does not exist yet. Writing..."
			}
		}
		NR%4==1{
			if ( ++readsProcessed%splitPoint == 0 ) {	# Multiple of splitPoint, increment files.
				blockCount++
				padBlockCount=sprintf("%0"zeroPad"d", blockCount)
				# Check if we are writing blocks or if we are not writing blocks then if read is 1, check for an alignment run.
				if (writeBlock || readNumber=="R1") {
					# Spawn alignment if the next block in the sequence exists for both reads.
					if (system(pBin"/check_blocks.sh "sampleID" "prefix" "readNumber" "curBlock" "blockCount" "alignArray" "sortArray) != 0) {
						print outHeader": CheckBlock failure! Aborting."
						exit 1
					}
				}
				
				if (writeBlock) {
					close (outStream)
					system("touch "outFile".done")
					print outHeader": Block "padBlock" finished at "readsProcessed" reads. Starting "padBlockCount
				} else {
					print outHeader": Block "padBlock" already written. Moving on to "padBlockCount
				}
				
				# Update current block number
				curBlock=blockCount
				padBlock=sprintf("%0"zeroPad"d", curBlock)
				outFile="blocks/"prefix"_"readNumber"_"padBlock".fastq.gz"
				
				if (system("[ -e "outFile".done ]") == 0) {
					# The new outfile.done already exists! Skip this block
					print outHeader": Skipping block "padBlock" as it is already completed!"
					writeBlock=0
				} else {
					writeBlock=1
				}
			}
			
			# Check if index sequence varies at too many positions.
			bestIndexChars=split(seqIndex,compIndex,"")
			indexChars=split($10,lineIndex,"")
			indexDelta=0;
			for (i=1; i<=indexChars; i++) {
				if (lineIndex[i] != compIndex[i]) {
					indexDelta++
				}
			}
			
			# If the index it too different, make a new file just for it outside the current sequence.
			if (indexDelta > maxDelta) {
				prefix=$2"_"$3"_"$4"_"$5"_"$10
			} else {
				prefix=$2"_"$3"_"$4"_"$5"_"seqIndex
			}
			
			# If read-group isn not within bounds, add to list.
			if (length(old_prefix) == 0 || prefix != old_prefix ) {
				print prefix > "blocks/"readNumber"_ReadGroup.txt"
			}
			old_prefix=prefix
			
			outFile="blocks/"prefix"_"readNumber"_"padBlock".fastq.gz"
			
			outStream=compCmd" > "outFile
		}
		
		{
			if (writeBlock) print | outStream
		}
		
		END{
			if (writeBlock) {
				close (outStream)
				system("touch "outFile".done")
				print outHeader": Block "blockCount" with "readsProcessed" read."
			} else {
				print outHeader": Block "blockCount" already exists with "readsProcessed" reads."
			}
			
			system(pBin"/check_blocks.sh "sampleID" "prefix" "readNumber" "curBlock" "curBlock" "alignArray" "sortArray)
		}
	'; then
		exit 1
	fi
}

[ "$(which pigz)" != "" ] && ZIP_CMD="pigz" || ZIP_CMD="gzip"
[ "${INPUT##*.}" != "gz" ] && CAT_CMD=cat || CAT_CMD="${ZIP_CMD} -cd"

echo "$HEADER: Zip command: ${ZIP_CMD}"
echo "$HEADER: Cat command: ${CAT_CMD}"

bestIndex=$(getBestIndex)

echo "$HEADER: Best index is [${bestIndex}]"

mkdir -p blocks

if ! splitByReadGroupAndCompress; then
	cmdFailed
	exit $EXIT_PR
fi

#rm ${INPUT} && echo "$HEADER: Purged input file!"

touch ${SAMPLE}_R${READNUM}_split.done
