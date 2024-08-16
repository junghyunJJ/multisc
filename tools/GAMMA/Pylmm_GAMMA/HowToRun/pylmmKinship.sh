#!/bin/tcsh
#$ -N esitmateK
#$ -cwd
#$ -o /u/home/j/jwjjoo/project/GAMMA/microbiota/Metisim//clusterJob/estimateK.$JOB_ID.out
#$ -j y
#$ -l h_rt=23:59:59
module load python
python /u/home/j/jwjjoo/project/program/pylmm/pylmmKinship.py --tfile=/u/home/j/jwjjoo/project/plink /u/home/j/jwjjoo/project/plink.K
python /u/home/j/jwjjoo/project/program/pylmm/pylmmKinship.py --emmaSNP=/u/home/j/jwjjoo/project/X.txt(snp by indi) --emmaNumSNPs=100 /u/home/j/jwjjoo/project/plink.K
