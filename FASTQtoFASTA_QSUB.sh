#!/bin/bash
#makeFLAShoutFile="$3.extendedFrags.fastq"
#FASTA="$3.fasta"
command="python /srv/gsfs0/projects/mignot/TCR_PIPELINE/FASTQtoFASTAv3.py /srv/gsfs0/projects/mignot/TCR_PIPELINE/$1 /srv/gsfs0/projects/mignot/TCR_PIPELINE/$2"
touch FASTQ_FASTA.sh
chmod 755 FASTQ_FASTA.sh
echo \#\!/bin/bash >FASTQ_FASTA.sh
echo \#$ -N FASTQ_to_FASTA >>FASTQ_FASTA.sh
echo \#$ -l h_vmem=80G >>FASTQ_FASTA.sh
echo \#$ -l h_rt=6:00:00 >>FASTQ_FASTA.sh
echo \#$ -hold_jid FLASH_TCR >>FASTQ_FASTA.sh
#echo \#$ -pe shm 8 >>FLASH_TEST.sh
echo "module load python/2.7.9" >>FASTQ_FASTA.sh
echo $command >> FASTQ_FASTA.sh
qsub -V -cwd FASTQ_FASTA.sh


