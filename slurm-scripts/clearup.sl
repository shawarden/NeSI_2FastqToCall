#!/bin/bash
#SBATCH --account		uoo00032
#SBATCH --job-name		ClearUp
#SBATCH --time			359
#SBATCH --mem-per-cpu	512
#SBATCH --cpus-per-task	6
#SBATCH --error			slurm/CU_%j.out
#SBATCH --output		slurm/CU_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

IDN=${1}

HEADER="CU"

echo "$HEADER: Saving work tree."

# Get remaining data store size
du -sh

# Purge excess files.
# This must be run AFTER the upload completes otherwise...
find . -type f -regex '.*\(bam\|bai\|vcf\|gz\|tbi\)$' -exec sh -c 'echo "$HEADER: Purging file $(realpath {})"; rm {}' \;

# Get data-less store size.
du -sh

# Roll up entire directory structure and .done files.
tar cf - ${PWD} | ${ZIP_CMD} -9 > ../${IDN}.tar.gz

# Upload workflow structure.
if ! . ${SLSBIN}/transfer.sl ${IDN} ../${IDN}.tar.gz; then
	echo "$HEADER: Transfer process failed!"
	exit $EXIT_TF
fi

# Purge run.
# This really should fail as the --error and --output files are in this folder.
rm --interactive=never -fr *
rm ../${IDN}.tar.gz;
