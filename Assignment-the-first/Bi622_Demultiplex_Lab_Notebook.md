Author: Evelyn Wong
Date Created: 2023-07-26
Last Updated: 2023-07-27

#2023-07-26
Performing initial data exploration. 
Leslie showed the intial number of lines in the files (over 2 billion lines) which took quite a while. 
We discussed in-class about what files are which:
Both R1 and R4 are the sequence reads, whereas R2 and R3 are the indexes. 

Because the files are HUGE, Leslie advised us that we should not unzip these files nor should we loop through the files more than once. 

1. Initial exploration! 
--------------------
Do not want to run on login node of cours (otherwise admin would be very angry)

Set up compute node.
srun --account=bgmp --partition=compute --time=1:00:00 --pty bash

Move into directory with the data.
cd /projects/bgmp/shared/2017_sequencing/ 

ii. Confirming number length of reads in each file:

(For the index files)
zcat 1294_S1_L008_R2_001.fastq.gz | head -n 2 | grep -v ^@ | tr -d "\n" | wc -c
Output: 8
zcat 1294_S1_L008_R3_001.fastq.gz | head -n 2 | grep -v ^@ | tr -d "\n" | wc -c
Output: 8
zcat 1294_S1_L008_R1_001.fastq.gz | head -n 2 | grep -v ^@ | tr -d "\n" | wc -c
Output: 101
zcat 1294_S1_L008_R4_001.fastq.gz | head -n 2 | grep -v ^@ | tr -d "\n" | wc -c
Output: 101

iii. Determine the phred encoding:

(Looking at the first 8 lines)
zcat 1294_S1_L008_R1_001.fastq.gz | sed -n 4~4p | grep -v ^@ | head -8
Output:
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ
A#AAFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJAJJJJJJJJJJJJJJFJJJJJFFFFJJJJJJJJJJJJJJJJJJ77F
A#AFFFJFJJFJJJJFJJJJJJJJAJJFJJJJJFJFJ7<FAFJJFJFJJFJFJJJFJAAJJJFJJJJJJJJJJJJJJJAJJJFAJJJJJFFJJJAJJJ<F-
A#<AAFJFJJJJFJJFJJ7JFJJJFJFAJJ<FF<<JJ<JJ<F<JJFAJJFFFJJJJJJA--77FJ--<<-AA<<AFJJJJJJFJJJFFFJ-<7--7-FFFA
A#A-FJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJFJFJFJJJJJJJJJJFFJJJJJJFJJJJJJJJJJJJJJJJJJFJ<FJJAFJJJF<J7FJJJF
-#AAFAAFJJFJJJJJJJJJJJJJF77<JJ<JFJAJFJFFA-F-A7JAJJJA<FAFFJJJFAFJF<FJJ-FAJFJJJJ<AJFJFFJJJJA-F-<AAA<FJJ
A#AFF<FJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJ<JAJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJFJJJJJJ-AJJJAJJFJJJFJF
A#AF<FFJJJFFJJJJJJJJJJJJJJJFJJJJJJJJJFFJJJJJJJJJJJJJJJJJAJAJJJJJFJJFJJJJJJJJJJJJJJ<JFJJJFJJFJFJJJJAJJ

Notice the '#' found only in Phred + 33. 

2. Generate a per base distribution of quality scores for read1, read2, index1, index2. For this, I'm reusing PS4 but without the functions. 

Create the list of 101 for the reads; create list of 8 for the indexes. 
Using argsparse to pass the length in my code. 
Created a slurm script
Keeps exiting with status 2 and/or 1. 
Figured it out: I forgot to not include spaces in my script for the file directories.
Fixed:
#dir1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
#dir2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
#dir3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
dir4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz

#/usr/bin/time -v /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-first/Bi622_Pt1_Qscore_Dist.py -f $dir1 -l 101 -o /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-first/Read1
#/usr/bin/time -v /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-first/Bi622_Pt1_Qscore_Dist.py -f $dir2 -l 8 -o /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-first/Index1
#/usr/bin/time -v /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-first/Bi622_Pt1_Qscore_Dist.py -f $dir3 -l 8 -o /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-first/Index2
/usr/bin/time -v /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-first/Bi622_Pt1_Qscore_Dist.py -f $dir4 -l 101 -o /projects/bgmp/evew/bioinfo/Bi622/Demultiplex/Assignment-the-first/Read2

While running each one, I'm commenting them out. Ran it on the small index files first. Takes around 18 minutes each. 
Ran both of the big read files, takes around 1:30:00 to 2:00:00. 

***I just read the assignment. Also wants the number of N's too? I'll have to remember to add that in next time...***

Have to make the unit tests! I forgot. Taking the head of the first few files records each for the "test" file. How should I add the headers for the unknown, matched, and unmatched...? Little confused. Are we doing this manually?
Asked Ross how he did his. He added his manually. I'll try to do the same here. 
So 6 files. Matched has 2 each. Unknown has 1. Unmatched has 1. 

Oh. Gitignore doesn't look at the fastq files. Have to rename my unit tests...whoops.


