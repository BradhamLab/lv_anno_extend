#! /bin/bash -l
#$ -N add_annos
#$ -M dyh0110@bu.edu
#$ -m eas
#$ -P bradham

source activate alignment;
snakemake --cluster 'qsub -P {cluster.project} ' --jobs 4 --latency-wait 30
