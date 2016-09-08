#!/bin/bash
#SBATCH --job-name		GenderDetermination
#SBATCH --time			0-00:10:00
#SBATCH --mem-per-cpu	512
#SBATCH --cpus-per-task	1
#SBATCH --error			slurm/GD_%j.out
#SBATCH --output		slurm/GD_%j.out

source /projects/uoo00032/Resources/bin/NeSI_2FastqToCall/baserefs.sh

IDN=${1}
PLATFORM=${2}
OUTPUT=coverage.sh

echo "#This file contains gender definitions for individual ${IDN}." | tee ${OUTPUT}
echo "#GD: ${PLATFORM}" | tee -a ${OUTPUT}
echo "" | tee -a ${OUTPUT}

source ${PLATFORMS}/${PLATFORM}.sh

   aCount=0
aCoverage=0
xCoverage=0
yCoverage=0

ren='^[0-9]+$'

printf "#%-10s %8s %s\n" "Chromosome" "Coverage" "Count" | tee -a ${OUTPUT}

for contig in ${CONTIGARRAY[@]}; do
	INPUT=depth/${contig}.sample_summary
	if [ ! -e ${INPUT} ]; then
		# Oh crappola!
		echo "#${contig} file ${INPUT} doesn't exist!" | tee -a ${OUTPUT}
		echo "exit 1" | tee -a ${OUTPUT}
		exit $EXIT_IO
	fi
	
	contigCount=$(awk 'NR==2{print $2}' ${INPUT})
	contigCover=$(awk 'NR==2{print $3}' ${INPUT})
	printf "#%-10s %8s %s\n" ${contig} ${contigCover} ${contigCount} | tee -a ${OUTPUT}
	
	if [[ $contig =~ $ren ]]; then	# Contig is a number
		aCount=$(($aCount + $contigCount))
		aCoverage=$(echo "($aCoverage + ($contigCount * $contigCover))" | bc)
	elif [[ $contig == X ]]; then
		xCoverage=${contigCover}
	elif [[ $contig == Y ]]; then
		yCoverage=${contigCover}
	fi
done

echo "" | tee -a ${OUTPUT}

aCoverage=$(echo "scale=2;$aCoverage/$aCount" | bc)

xaRatio=$(echo "scale=3; $xCoverage/$aCoverage" | bc | sed 's/^\./0./')
yaRatio=$(echo "scale=3; $yCoverage/$aCoverage" | bc | sed 's/^\./0./')
xCount=$(echo "scale=3; ($xCoverage/$aCoverage)/$XRat" | bc | sed 's/^\./0./')
yCount=$(echo "scale=3; ($yCoverage/$aCoverage)/$YRat" | bc | sed 's/^\./0./')

printf "%-20s %6s\n" "#Autosomal coverage:" ${aCoverage} | tee -a ${OUTPUT}
printf "%-20s %6s\n" "#True-X coverage:" ${xCoverage} | tee -a ${OUTPUT}
printf "%-20s %6s\n" "#True-Y coverage:" ${yCoverage} | tee -a ${OUTPUT}
printf "%-20s %6s\n" "#X:A ratio:" ${xaRatio} | tee -a ${OUTPUT}
printf "%-20s %6s\n" "#Y:A ratio:" ${yaRatio} | tee -a ${OUTPUT}
printf "%-20s %6s\n" "#Fractional X:" ${xCount} | tee -a ${OUTPUT}
printf "%-20s %6s\n" "#Fractional Y:" ${yCount} | tee -a ${OUTPUT}

Xmin=$(echo "scale=3;$xCount - $XVar" | bc | sed 's/^\./0./')
Xmax=$(echo "scale=3;$xCount + $XVar" | bc | sed 's/^\./0./')
Ymin=$(echo "scale=3;$yCount - $YVar" | bc | sed 's/^\./0./')
Ymax=$(echo "scale=3;$yCount + $YVar" | bc | sed 's/^\./0./')

xChromes=0
yChromes=0

# Count X and Y chromosomes that fall within boundries from whole numbers between 1 and 4: XXXXYYYY at most.
for i in `seq 4`
do
	XinRange=$(echo "$i <= $Xmax && $i >= $Xmin" | bc)
	YinRange=$(echo "$i <= $Ymax && $i >= $Ymin" | bc)
	
	if [ $XinRange -eq 1 ]
	then
		printf "%-20s %6s %s\n" "#X (boundry):" $i "($Xmin < $xCount < $Xmax)" | tee -a ${OUTPUT}
		xChromes=$i
	fi
	
	if [ $YinRange -eq 1 ]
	then
		printf "%-20s %6s %s\n" "#Y (boundry):" $i "($Ymin < $yCount < $Ymax)" | tee -a ${OUTPUT}
		yChromes=$i
	fi
done


# Build chromosome string.
if [ $xChromes -gt 0 ]
then
	# There are X chromosomes within the defined boundries.
	# Write a line of that many Xs.
	sexchromosomes=$(for a in `seq ${xChromes}`; do echo -n X; done)
elif [ $(echo "scale=3;$xCount > (1.0 - $XVar)" | bc) -eq 1 ]
then
	# There are no X chromosomes within the boundries.
	# Frational portions of X are greater than ONE.
	# Append these chromosomes with an E mark!
	sexchromosomes=E$(for a in `seq ${xCount}`; do echo -n X; done)
else
	# There are no X chromosomes within the boundries.
	# Frations of X found are below the lowest possible boundry.
	# Set the number of X chromoromes to ZERO.
	sexchromosomes="0"
fi

if [ $yChromes -gt 0 ]
then
	# There are Y chromosomes within the defined boundries.
	# Write a line of that many Ys
	sexchromosomes=${sexchromosomes}$(for a in `seq ${yChromes}`; do echo -n Y; done)
elif [ $(echo "scale=3;$yCount > (1.0 - $YVar)" | bc) -eq 1 ]
then
	# There are no Y chromosomes within the boundries.
	# Fraction portions of Y are greater than ONE.
	# Append these chromomes with an E mark.
	sexchromosomes=${sexchromosomes}E$(for a in `seq ${yCount}`; do echo -n Y; done)
elif [ $xChromes -eq 1 ]
then
	# There are no Y chromosomes within the boundries.
	# There is ONE X chromosome.
	# Fractional portions of Y are below lowest possible boundry.
	# Set the number of Y chromosomes to ZERO.
	sexchromosomes=${sexchromosomes}0
fi

# Decide overall gender
if [[ $xChromes -eq 0 ]]
then
	# Could not find any X chromosomes so FAIL!
	calculatedgender="Unknown"
else
	# There is at least ONE X chromosome
	if [[ $yChromes -eq 0 ]] && [[ $(echo "scale=3;$yCount < (1.0 - $YVar)" | bc) -eq 1 ]]
	then
		# There are no Y chromosomes within the boundries.
		# There are no fractional Y chromosome portions greater than the lowest possible broundry.
		calculatedgender="Female"
	else
		# There are at least 1 full Y chromosome present, even if it falls outside the boundries.
		calculatedgender="Male"
	fi
fi

printf "%-20s %6s\n" "#SexChr:" $sexchromosomes | tee -a ${OUTPUT}
printf "%-20s %6s\n" "#Gender:" $calculatedgender | tee -a ${OUTPUT}

if [ $xChromes -eq 0 ]; then
	echo "#Unknown gender. Please check capture platform or sample contamination!" | tee -a ${OUTPUT}
	echo "exit 1" >> ${OUTPUT}
	exit $EXIT_PR
fi

echo "" | tee -a ${OUTPUT}

echo "X_CHROMOSOMES=$xChromes" | tee -a ${OUTPUT}
echo "Y_CHROMOSOMES=$yChromes" | tee -a ${OUTPUT}

touch ${OUTPUT}.done

if ! ${SLSBIN}/transfer.sl ${IDN} coverage.sh ${IDN}.coverage.sh; then
	echo "$HEADER: Transfer failed!"
	exit $EXIT_TF
fi
