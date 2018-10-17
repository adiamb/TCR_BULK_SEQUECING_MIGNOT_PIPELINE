#!/bin/bash
command="rm $1.notCombined_1.fastq $1.notCombined_2.fastq $1.extendedFrags.fastq $1.hist $1.hist.innie $1.hist.outie $1.histogram $1.histogram.innie $1.histogram.outie "
touch FLASH_FILE_RM.sh
chmod 755 FLASH_FILE_RM.sh
echo \#\!/bin/bash >FLASH_FILE_RM.sh
echo \#$ -N rem_FLASH_files >>FLASH_FILE_RM.sh
echo \#$ -l h_rt=6:00:00 >>FLASH_FILE_RM.sh
echo \#$ -l h_vmem=10G >>FLASH_FILE_RM.sh
echo \#$ -hold_jid clumping_rcl_TCR >>FLASH_FILE_RM.sh
echo $command >> FLASH_FILE_RM.sh
qsub -V -cwd FLASH_FILE_RM.sh
