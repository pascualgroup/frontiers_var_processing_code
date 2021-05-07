# frontiers paper road map

Qixin He

This document explains the simulation code, experimental design and code to analyze and generate the figures.

## Simulation code: 
`varmodel2` version e0c619b
(https://github.com/pascualgroup/varmodel2/tree/e0c619b975ac5fcc6b709c2204e7395f7197a636)
There are newer commits for Frederic's project. Do not use those.

## Scheduling code
all the scheduling codes, biting rates pre- and post IRS, biting rates during IRS are stored in the folder `IRStests` in `varmodel2`

1. Scenario settings in `AugNew_param.csv`. Each column depicts a parameter that would vary per run

2. Input file template file: `varMigTestTemplate.py`. Parameters in brackets will be replaced by actual values per run.

3. python script to write the actual input file according to the template and scenario No.: `writeGIParameters.py`

This python script outputs at least 8 input files for each scenario:
s0-s3: specific selection
s4-s7: generalized immunity

s0: preIRS, specific selection
s1,s2,s3: IRS=2,5,and 10years
s4: preIRS, generalized immunity
s5-s7: IRS=2,5,and 10years

`-r` option determines how many replicates are ran per scenario. In this study, we put `-r=50`

To run the script for scenario 13 with the run name `AugNew`, type:
`python writeGIParameters.py -p AugNew_param.csv -i varMigTestTemplate.py -n 13 -r 50 -x AugNew`

## Code to analyze the sqlite output

`writeSummaryTable.r`, which requires `networkCalculations.r`
This code takes the input sqlite, write out a table with basic summary statistics, genetic diversity measures, network property calculations, pts distributions, and doi distribution.

## sbatch file to submit to the cluster:

```
module load gcc/6.1
module load python/3.7.0
module load R

filePrefix=AugNew

#here define where the actual run will be stored
remoteDir=/scratch/midway2/heqixin/${filePrefix}_${SLURM_ARRAY_TASK_ID}
#copy original Rdiv code to the running folder
cp -r /home/heqixin/varmodel2 $remoteDir
```

generate all the input files for the scenario

```
cd $remoteDir
cd IRStests
python writeGIParameters.py -p ${filePrefix}_param.csv -i varMigTestTemplate.py -n $SLURM_ARRAY_TASK_ID -r 50 -x $filePrefix
cd ..
```

build the model, run preIRS first
```
./build.py -p IRStests/${filePrefix}_${SLURM_ARRAY_TASK_ID}_s0_input.py -d s0
#execute the run
./s0/bin/varMig
Rscript ../runInfo/writeSummaryTable.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ 0 0
```

run GI next
```
cd IRStests
python writeGIParameters.py -p ${filePrefix}_param.csv -i varMigTestTemplate.py -n $SLURM_ARRAY_TASK_ID -g -r 50 -x $filePrefix
cd ..

./build.py -p IRStests/${filePrefix}_${SLURM_ARRAY_TASK_ID}_s4_input.py -d s4
#execute the run
./s4/bin/varMig

Rscript ../runInfo/writeSummaryTable.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ 4 0
```

execute IRS runs

```
for j in {0..49}
do
	for i in {1..3}
	do
	./build.py -p IRStests/${filePrefix}_${SLURM_ARRAY_TASK_ID}_s${i}_r${j}_input.py -d s${i}_r${j}
	./s${i}_r${j}/bin/varMig
	Rscript ../runInfo/writeSummaryTable.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ ${i} ${j}
	done

	for i in {5..7}
	do
	./build.py -p IRStests/${filePrefix}_${SLURM_ARRAY_TASK_ID}_s${i}_r${j}_input.py -d s${i}_r${j}
	./s${i}_r${j}/bin/varMig
	Rscript ../runInfo/writeSummaryTable.r ./ ${filePrefix} ${SLURM_ARRAY_TASK_ID} ../results/ ${i} ${j}
	done
done
```

## post processing scripts and generate graphs

1. first combine all the results output to one file:
`cat *_summaryTable.txt > [filePrefix]_allSumTab.txt`
`cat *_ptsTable.txt >[filePrefix]_allPTS.txt` 
`cat *_netTable.txt >[filePrefix]_allNetCal.txt` 
`cat *_durTable.txt >[filePrefix]_allDurTab.txt` 

2. perform analysis and plot figures

see comments in `analysis_code_varMig.R` to follow the steps.





