#!/bin/bash
#SBATCH --account uoo00032
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sam.hawarden@otago.ac.nz
#SBATCH --mail-type=FAIL
#SBATCH --constraint=avx
#SBATCH --error=slurm/catbam_%j.out
#SBATCH --output=slurm/catbam_%j.out

source /projects/uoo00032/Resources/bin/baserefs.sh

INPUT="printreads/1/3507.bam"
OUTPUT="samtools_cat_3507.bam"
HEADER="printreads/1/3507.bam"

echo "SamView: Header ${HEADER} + ${INPUT} -> ${OUTPUT}"
date

for file in ${INPUT}; do
	if [ ! -e ${file} ]; then
		echo "SamView: \"${file}\" doesn't exist!"
		scriptFailed "SamView"
		exit 1
	fi
done

inputCount=$(echo ${INPUT} | wc -w)
contigCount=$(echo ${CONTIGS} | wc -w)

module load ${MOD_SAMTOOLS}

CMD="$(which srun) $(which samtools) view -H ${INPUT}"
echo "SamView: ${CMD}" | tee -a commands.txt

${CMD}
passed=$?

echo "SamView: ${INPUT} -> ${OUTPUT} ran for $(($SECONDS / 3600))h $((($SECONDS % 3600) / 60))m $(($SECONDS % 60))s"
date

if [ $passed -ne 0 ]
then
	echo "SamView: ${INPUT} -> ${OUTPUT} Failed!"
	scriptFailed "SamView"
	exit 1
fi

touch ${OUTPUT}.done

storeMetrics
