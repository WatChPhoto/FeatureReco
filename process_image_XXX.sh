#!/bin/bash
#SBATCH --account=rpp-blairt2k
#SBATCH --time=02:00:00
#SBATCH --mem=12000M
#SBATCH --cpus-per-task=11
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load gcc/7.3.0
module load nixpkgs/16.09
module load opencv/4.4.0
module load root/6.14.04


AAAA=`grep "^XXX " circleranges.txt | awk '{print $2}'`
BBBB=`grep "^XXX " circleranges.txt | awk '{print $3}'`

cd /scratch/blairt2k/process_images/BarrelFar
mkdir XXX
cp Config.txt XXX/.
cd XXX
sed -i "s/AAAA/${AAAA}/g" Config.txt
sed -i "s/BBBB/${BBBB}/g" Config.txt
echo "Starting to run"
/scratch/blairt2k/FeatureReco/FindBoltLocations /scratch/blairt2k/BarrelSurveyFar_TB3/BarrelSurveyFar_median/images/XXX.jpg
mv FindBoltLocation.root FindBoltLocationXXX.root
echo "Done execution"


