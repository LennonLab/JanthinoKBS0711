#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load fastqc
module load bioperl
module load python
module load cutadapt
module load khmer
module load velvet
module load VelvetOptimiser
PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
cd /N/dc2/projects/Lennon_Sequences/2014_Janthino

cutadapt -q 10 -a AGATCGGAAGAGC --minimum-length 20 -o ./Genom_Announ_Test/tmp.1.fastq -p ./Genom_Announ_Test/tmp.2.fastq ./KBS0711_wt_R1.fastq ./KBS0711_wt_R2.fastq
cutadapt -q 10 -a AGATCGGAAGAGC --minimum-length 20 -o ./Genom_Announ_Test/KBS0711_wt_R2_no_adapt.fastq -p ./Genom_Announ_Test/KBS0711_wt_R1_no_adapt.fastq ./Genom_Announ_Test/tmp.2.fastq ./Genom_Announ_Test/tmp.1.fastq

cd Genom_Announ_Test/
rm tmp.1.fastq tmp.2.fastq
interleave-reads.py ./KBS0711_wt_R1_no_adapt.fastq ./KBS0711_wt_R2_no_adapt.fastq -o ./KBS0711_wt_no_adapt.fastq      
fastq_quality_filter -Q33 -q 30 -p 50 -i ./KBS0711_wt_no_adapt.fastq > ./KBS0711_wt_no_adapt.trim.fastq
extract-paired-reads.py ./KBS0711_wt_no_adapt.trim.fastq
normalize-by-median.py -k 25 -C 25 -N 4 -x 2e9 -p ./KBS0711_wt_no_adapt.trim.fastq.pe
extract-paired-reads.py ./KBS0711_wt_no_adapt.trim.fastq.pe.keep
VelvetOptimiser.pl -s 31 -e 68 -f '-shortPaired -fastq ./KBS0711_wt_no_adapt.trim.fastq.pe.keep.pe' -t 8 -k 'n50*ncon' --p 'wt_janthin_genom_ann'
