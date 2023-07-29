#!/usr/bin/env python
# Author: Evelyn Wong
# Date Created: 2023-07-27
# Last Updated: 2023-07-28

import argparse
import gzip
import bioinfo
import matplotlib.pyplot as plt 

def get_args():
    '''This function gets arguments for FASTA file from user to be parsed'''
    parser = argparse.ArgumentParser(description="This function retrieves the FASTA file of interest to be parsed")
    parser.add_argument("-f", "--file_name", help = "Enter file name to be parsed", required = True)
    parser.add_argument("-l", "--read_length", type = int, help = "Enter read length", required = True)
    parser.add_argument("-o", "--output_file", help = "Enter output file name", required = True)
    args = parser.parse_args()
    return args

args = get_args()
file_name = args.file_name
sequence_length: int = args.read_length

# create a list with sequence length number of values inside
my_list: list = []
my_list = [0] * sequence_length
#my_list.extend([0] * sequence_length)
l_counter: int = 0 # counter to keep track of total num lines in file


# open the fastq file 
with gzip.open(file_name, "rt") as fh:
    for line in fh:
       #loop through the record
       l_counter += 1
       line = line.strip("\n")

       if l_counter % 4 == 0: 
            for num in range(len(line)):
                my_list[num] +=  bioinfo.convert_phred(line[num])
    
for base in range(len(my_list)):
    my_list[base] = my_list[base]/(l_counter/4)

plt.bar(range(len(my_list)), my_list)
plt.title('Mean Sums of Qscores per Base of Illumina Reads')
plt.xlabel('Q-score Base Positions')
plt.ylabel('Mean Q-score of Illumina Reads')
plt.savefig(f'{args.output_file}.png')