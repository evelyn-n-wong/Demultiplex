# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here: https://github.com/evelyn-n-wong/Demultiplex/blob/master/Assignment-the-first/Part_1/Bi622_Pt1_Qscore_Dist.py

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read1 | 101 | 33 |
| 1294_S1_L008_R2_001.fastq.gz | Index1 | 8 | 33 |
| 1294_S1_L008_R3_001.fastq.gz | Index2 | 8 | 33 |
| 1294_S1_L008_R4_001.fastq.gz | Read2 | 101 | 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. A good cut-off score for indices, because they are much smaller (only 8 nt), would be Q35 per base. If they have even one poor score listed, then they should be discarded or placed into the "unknown" bin as it would affect the entire sequence. For the reads, a cut-off score of ~30 for the mean position of all bases. This seems to be consistent with the graph as well.  
    3. Command used: ls -1 ./1294_S1_L008_R[2-3]_001.fastq.gz | while read file; do echo $file; zcat $file | sed -n 2~4p | grep "N" | wc -l; done
        ./1294_S1_L008_R2_001.fastq.gz
        3976613
        ./1294_S1_L008_R3_001.fastq.gz
        3328051
        Total: 7304664
    
## Part 2
1. Define the problem:
We want to look through a lane of sequencing generated from the 2017 BGMP cohort’s library preps and determine the level of index swapping and undetermined index-pairs, before and after quality filtering of index reads.
2. Describe output
Output: 24 matched indices for R1, 24 matched indices for R4, 2 files for unknown indices or low quality for R1 and R4, and 2 files for non-matching pairs for R1 and R4.
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode: 
https://github.com/evelyn-n-wong/Demultiplex/blob/master/Assignment-the-first/Part_2/Demultiplex_pseudocode.txt
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
