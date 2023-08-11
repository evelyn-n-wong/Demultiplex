#!/bin/bash
#Author: Evelyn Wong
#Date Created: 2023-08-07

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=compute               #REQUIRED: which partition to use
#SBATCH --mail-user=evew@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB

conda activate bgmp_py311 

dir1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
dir2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
dir3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
dir4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
indices=/projects/bgmp/shared/2017_sequencing/indexes.txt

#/usr/bin/time -v /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-third/Test_code.py -r1 $dir1 -r2 $dir2 -r3 $dir3 -r4 $dir4 -i $indices
/usr/bin/time -v /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-third/demultiplex.py -r1 $dir1 -r2 $dir2 -r3 $dir3 -r4 $dir4 -i $indices
