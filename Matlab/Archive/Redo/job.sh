#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=100:0:0
#$ -R y
#$ -pe smp 5

#$ -cwd

#$ -S /bin/bash
#$ -j y
#$ -N testMC

hostname
date

source /share/apps/source_files/matlab

date
matlab -nodisplay -nosplash -nodesktop -r "run('Simulation.m');exit;"
date
