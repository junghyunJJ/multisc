#!/bin/tcsh
#$ -N Project_1
#$ -cwd
#$ -o /u/home/j/jwjjoo/project/clusterJob/project.$JOB_ID.out
#$ -j y
#$ -l h_rt=23:59:59
module load python
python /u/home/j/jwjjoo/project/program/pylmm/pylmmGWAS.py -v --tfile=/u/home/j/jwjjoo/project/plink --kfile=/u/home/j/jwjjoo/project/plink.K --phenofile=/u/home/j/jwjjoo/project/project_1.pheno /u/home/j/jwjjoo/project/clusterJob/project_1.log
python /u/home/j/jwjjoo/project/program/pylmm/pylmmGWAS.py -v --emmaSNP=/u/home/j/jwjjoo/project/X.txt(snp by indi) --kfile=/u/home/j/jwjjoo/project/plink.K --emmaPHENO=/u/home/j/jwjjoo/project//pheno_t.txt(one row file) /u/home/j/jwjjoo/project/clusterJob/test.log
python /u/home/j/jwjjoo/project/program/pylmm/pylmmGWAS_multiPhHeri.py -v --emmaSNP=/u/home/j/jwjjoo/project/X.txt --kfile=/u/home/j/jwjjoo/project/plink.K --emmaPHENO=/u/home/j/jwjjoo/project/Y.txt(pheno by indi) /u/home/j/jwjjoo/project/clusterJob/test.log
