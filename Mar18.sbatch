#!/bin/bash
#SBATCH --job-name=Mar18
#SBATCH --time=36:00:00
#SBATCH --output=/scratch/midway2/heqixin/runInfo/Mar18_%A_%a.out
#SBATCH --error=/scratch/midway2/heqixin/runInfo/Mar18_%A_%a.err
#SBATCH --array=14,27
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --partition=broadwl
#SBATCH --account=pi-pascualmm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=heqixin@uchicago.edu

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
module load gcc/6.1
module load python/cpython-3.7.0
module load R

filePrefix=LongEnd

#here define where the actual run will be stored
remoteDir=/scratch/midway2/heqixin/${filePrefix}_${SLURM_ARRAY_TASK_ID}
#copy original Rdiv code to the running folder
cp -r /home/heqixin/varmodel2 $remoteDir

cd $remoteDir
cd IRStests
python writeGIParameters.py -p ${filePrefix}_param.csv -i {filePrefix}_Template.py -n $SLURM_ARRAY_TASK_ID -r 5 -x $filePrefix
cd ..
#build the model, run preIRS first
./build.py -p IRStests/${filePrefix}_${SLURM_ARRAY_TASK_ID}_s0_input.py -d s0
#execute the run
./s0/bin/varMig
Rscript writeSummaryTableFast.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ 0 0
#Rscript writeSummaryTable.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ 0 0

#run GI next
cd IRStests
python writeGIParameters.py -p -p ${filePrefix}_param.csv -i {filePrefix}_Template.py -n $SLURM_ARRAY_TASK_ID -g -r 5 -x $filePrefix
cd ..

./build.py -p IRStests/${filePrefix}_${SLURM_ARRAY_TASK_ID}_s4_input.py -d s4
#execute the run
./s4/bin/varMig

Rscript writeSummaryTableFast.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ 4 0
#Rscript writeSummaryTable.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ 4 0



#execute IRS runs
for j in {0..4}
do
	for i in {3..3}
	do
	./build.py -p IRStests/${filePrefix}_${SLURM_ARRAY_TASK_ID}_s${i}_r${j}_input.py -d s${i}_r${j}
	./s${i}_r${j}/bin/varMig
	Rscript writeSummaryTableFast.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ ${i} ${j}
	Rscript ../runInfo/writeDuration.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ ${i} ${j}
	
	done

	for i in {7..7}
	do
	./build.py -p IRStests/${filePrefix}_${SLURM_ARRAY_TASK_ID}_s${i}_r${j}_input.py -d s${i}_r${j}
	./s${i}_r${j}/bin/varMig
	Rscript writeSummaryTableFast.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ ${i} ${j}
	Rscript ../runInfo/writeDuration.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ ${i} ${j}
	done
done
	

