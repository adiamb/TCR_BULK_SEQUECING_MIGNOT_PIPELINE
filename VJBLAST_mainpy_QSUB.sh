#!/bin/bash
command="python /srv/gsfs0/projects/mignot/TCR_PIPELINE/TCR_VJ_ID_BLAST_Analysis_v5_5mer_AA.py
/srv/gsfs0/projects/mignot/TCR_PIPELINE/$1_chunk\$SGE_TASK_ID.fasta 
/srv/gsfs0/projects/mignot/TCR_PIPELINE/TCR_V_V2_complement_trimmed_DB
/srv/gsfs0/projects/mignot/TCR_PIPELINE/TCR_J_V2_DB
/srv/gsfs0/projects/mignot/TCR_PIPELINE/TCR_sample5_DB
/srv/gsfs0/projects/mignot/TCR_PIPELINE/TCR_Identifier_updated.csv
/srv/gsfs0/projects/mignot/TCR_PIPELINE/$1_chunk\$SGE_TASK_ID.rcl
\$SGE_TASK_ID
/srv/gsfs0/projects/mignot/TCR_PIPELINE/TCR_C_DB
$1_TEST"
touch VJBLAST_Main.sh
chmod 755 VJBLAST_Main.sh
echo \#\!/bin/bash >VJBLAST_Main.sh
echo \#$ -N VJ_BLAST_TCR >>VJBLAST_Main.sh
echo \#$ -l h_vmem=8G >>VJBLAST_Main.sh
echo \#$ -l h_rt=6:00:00 >>VJBLAST_Main.sh
echo \#$ -pe shm 2 >>VJBLAST_Main.sh
echo \#$ -t 1-$2 >>VJBLAST_Main.sh
echo \#$ -hold_jid splitFASTATCRS >>VJBLAST_Main.sh
echo "module load blast/2.3.0+" >>VJBLAST_Main.sh
echo "module load python/2.7.9" >>VJBLAST_Main.sh
echo $command >> VJBLAST_Main.sh 
qsub -V -cwd VJBLAST_Main.sh

./VJBLAST_mainpy_6merBC_QSUB.sh