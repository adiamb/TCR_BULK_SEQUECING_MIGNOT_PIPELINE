#!/bin/bash
command="flash /srv/gsfs0/projects/mignot/TCR_PIPELINE/$1 /srv/gsfs0/projects/mignot/TCR_PIPELINE/$2 -m 5 -M 300 -O -o $3 -t 16"
touch FLASH_TEST.sh
chmod 755 FLASH_TEST.sh
echo \#\!/bin/bash >FLASH_TEST.sh
echo \#$ -N FLASH_TCR >>FLASH_TEST.sh
echo \#$ -l h_vmem=8G >>FLASH_TEST.sh
echo \#$ -l h_rt=6:00:00 >>FLASH_TEST.sh
echo \#$ -w e >>FLASH_TEST.sh
echo \#$ -pe shm 16 >>FLASH_TEST.sh
echo "module load flash/1.2.11" >>FLASH_TEST.sh
echo $command >> FLASH_TEST.sh
qsub -V -cwd FLASH_TEST.sh 

