#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=01:30:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/RI_%j.out
#SBATCH --output=slurm/RI_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

INPUT=${1}
OUTPUT=$([ "${2}" == "" ] && echo -ne "${INPUT%.bam}.bai" || echo -ne "${2}")
IDN=$(echo $SLURM_JOB_NAME | cut -d'_' -f2)

HEADER=RI

echo "$HEADER: ${INPUT} -> ${OUTPUT}"
date

# Make sure input and target folders exists and that output file does not!
if ! inFile;  then exit 1; fi
if ! outDirs; then exit 1; fi
if ! outFile; then exit 1; fi

CMD="$(which srun) ${SAMTOOLS} index ${INPUT} ${JOB_TEMP_DIR}/${OUTPUT}"
echo "$HEADER: ${CMD}" | tee -a commands.txt

if ! ${CMD}; then
	cmdFailed
	exit 1
fi

# Move output to final location
if ! finalOut; then exit 1; fi

touch ${OUTPUT}.done

sbatch -J TRI_${IDN} ${SLSBIN}/transfer.sl ${IDN} ${OUTPUT}

storeMetrics


