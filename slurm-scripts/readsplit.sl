#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/readsplit_%j.out
#SBATCH --output=slurm/readsplit_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

  SAMPLE=${1}
 READNUM=${2}
READFILE=${3}

echo "ReadSplit: ${READNUM} ${READFILE}"
date

if [ ! -e ${READFILE} ]; then
	echo "ReadSplit: Can't find ${READFILE}"
#	scriptFailed "ReadSplit"
	exit 1
fi

module load ${MOD_ZLIB}

# Scan the FastQ file for index sequence.
# Returns the index with the highest count within FASTQ_MAXSCAN lines from top.
function getBestIndex {
	${CAT_CMD} ${READFILE} | head -${FASTQ_MAXSCAN} | awk \
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
	${CAT_CMD} ${READFILE} | awk -F'[@:]' \
		-v sampleID="${SAMPLE}" \
		-v readNumber="${READNUM}" \
		-v seqIndex="${bestIndex}" \
		-v maxDelta="${FASTQ_MAXDIFF}" \
		-v splitPoint="${FASTQ_MAXREAD}" \
		-v compCmd="${ZIP_CMD}" \
		-v pBin="${PBIN}" \
		'
		BEGIN{
			blockCount=0
			curBlock=sprintf("%05d", blockCount)
			if (system("[ -e blocks/*_"readNumber"_"curBlock".fastq.gz.done ]") == 0) {
				# a blocks/..._Rn_00000.fastq.gz.done file exists so skip this block
				writeBlock=0
				print "ReadSplit: Skipping block "curBlock" as it is already completed!"
			} else {
				writeBlock=1
				print "ReadSplit: Block "curBlock" does not exist yet. Writing..."
			}
		}
		NR%4==1{
			if ( ++readsProcessed%splitPoint == 0 ) {	# Multiple of splitPoint, increment files.
				blockCount++
				if (writeBlock) {
					close (outStream)
					system("touch "outFile".done")
					print "ReadSplit: Block "curBlock" with "readsProcessed" reads. Now writing to "sprintf("%05d", blockCount)
				} else {
					print "ReadSplit: Block "curBlock" already written. Moving on to "sprintf("%05d", blockCount)
				}
				
				# Spawn alignment if the next block in the sequence exists for both reads.
				system("sleep 1s; "pBin"/check_blocks.sh "sampleID" "prefix" "readNumber" "curBlock" "sprintf("%05d", blockCount))
				
				# Update current block number
				curBlock=sprintf("%05d", blockCount)
				outFile="blocks/"prefix"_"readNumber"_"curBlock".fastq.gz"
				
				if (system("[ -e "outFile".done ]") == 0) {
					# The new outfile.done already exists! Skip this block
					print "ReadSplit: Skipping block "curBlock" as it is already completed!"
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
			
			outFile="blocks/"prefix"_"readNumber"_"curBlock".fastq.gz"
			
			outStream=compCmd" > "outFile
			
			if (writeBlock) print | outStream
		}
		
		NR%4==2{if (writeBlock) print | outStream}
		
		NR%4==3{if (writeBlock) print | outStream}
		
		NR%4==0{if (writeBlock) print | outStream}
		
		END{
			if (writeBlock) {
				close (outStream)
				system("touch "outFile".done")
				print "ReadSplit: Block "blockCount" with "readsProcessed" read."
			} else {
				print "ReadSplit: Block "blockCount" already exists with "readsProcessed" reads."
			}
			
			system(pBin"/check_blocks.sh "sampleID" "prefix" "readNumber" "curBlock" "curBlock)
		}'
}

[ "$(which pigz)" != "" ] && ZIP_CMD="pigz" || ZIP_CMD="gzip"
[ "${READFILE##*.}" != "gz" ] && CAT_CMD=cat || CAT_CMD="${ZIP_CMD} -cd"

echo "ReadSplit: Zip command: ${ZIP_CMD}"
echo "ReadSplit: Cat command: ${CAT_CMD}"

bestIndex=$(getBestIndex)

echo "ReadSplit: Best index is [${bestIndex}]"

mkdir -p blocks

splitByReadGroupAndCompress
passed=$?

echo "ReadSplit: Splitting ${READFILE} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]; then
	echo "ReadSplit: Failed to split ${READFILE}"
#	scriptFailed "ReadSplit"
	exit 1
fi

#if [ $(cat ${READNUM}_ReadGroup.txt | wc -w) -gt 1 ]
#then
#	for prefix in $(cat ${READNUM}_ReadGroup.txt)
#	do
#		echo "ReadSplit: Compressing ReadGroup ${prefix}: $(sbatch -J ReadMultiple_${prefix} s0_split_multiple.sl ${prefix}_${READNUM} | awk '{print $4}')"
#	done;
#	exit 1
#fi

rm ${READFILE}

touch ${SAMPLE}_${READNUM}_split.done

storeMetrics run
