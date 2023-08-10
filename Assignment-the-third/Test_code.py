#!/usr/bin/env python
# Author: Evelyn Wong
# Date Created: 2023-08-02
# Last Updated: 2023-08-03
# importing modules
import argparse
import gzip
import matplotlib.pyplot as plt
import itertools
import bioinfo

def get_args():
    '''This function gets arguments for FASTA file from user to be parsed'''
    parser = argparse.ArgumentParser(description="This function demultiplexes Illumina paired-end reads by index reads and outputs a summarized report for the user")
    parser.add_argument("-r1", "--read_one", help = "Enter R1 FASTQ file", required = True)
    parser.add_argument("-r2", "--index_one", help = "Enter R2 FASTQ file", required = True)
    parser.add_argument("-r3", "--index_two", help = "Enter R3 FASTQ file", required = True)
    parser.add_argument("-r4", "--read_two", help = "Enter R4 FASTQ file", required = True)
    parser.add_argument("-i", "--index_file", help = "Enter known indexes file", required = True)
    args = parser.parse_args()
    return args

args = get_args()

# Initalize variables for files
read1 = args.read_one # R1
index1 = args.index_one # R2
index2 = args.index_two # R3
read2 = args.read_two # R4
index_file = args.index_file

index_seq = set() # store all known indices
#known_indices: dict  = {} # dictionary to store all known indices from passed file
matched_in: dict = {} # dictionary to store all 24 matched, with value as integer with that matched
hopped_in: dict = {} # dictionary to store unmatched 
unknown_in: dict = {} # dictionary to store unknown and low-quality indices
score: int = 0 # qscore condition of each base 

# Initialize counters

num_matched: int = 0 # number of matched indices
num_unknown: int = 0 # number of unknown indices
num_low_qual: int = 0 # number of low quality indices
num_unmatched: int = 0 # number of unmatched indices
l_counter: int = 0 # number of lines
r_counter: int = 0 # number of records

with open(index_file, "r") as index_fh:
    # Opening the known index file and retrieving the index sequence for reference. 
    numLines: int = 0
    for line in index_fh:
        if numLines!= 0: # skip header
            line = line.strip("\n")
            index_columns = line.split()[4] # split the index columns and get only the sequence
            index_seq.add(index_columns) # add it to set
        numLines += 1
print(index_seq)

        # want viable index-pair combinations later on while writing 

def rev_complement(sequence: str):
    '''This function takes in a sequence and returns the reverse complement of the sequence as a string'''
    rc_dict: dict = {"A":"T", "T":"A", "G":"C", "C":"G","N":"N"}
    rev_comp: str = ""
    for char in sequence:
        rev_comp += rc_dict[char]
    rev_comp = rev_comp[::-1] #return reversed sequence 
    return rev_comp

def qual_check(qscore_line: str):
    '''This function takes in a quality score sequence, checks whether each base is not low-quality (less than 30 or average of whole line),
    and returns False if it fails condition'''
    qual_test: bool = True
    #q_average = bioinfo.qual_score(qscore_line) #get the average quality score for the entire q-score line
    for q in qscore_line: #go through base-by-base of the index
        score = bioinfo.convert_phred(q)
        #if score < q_average: #if the base of index is lower than average quality score of entire line
        if score < 30: 
            qual_test = False
            break
        return qual_test

def append_to_header(old_header: str, i1, i2):
    '''This function adds index-pair to header line for R1 and R4'''
    # Returns updated header consisting of old header line as a string, index 1, rev_comp of index 2
    new_header = old_header + ":" + i1 + "-" + i2
    index_pair = i1 + "-" + i2
    return new_header, index_pair

def create_record(file):
    '''This function reads in a file and gets the headers, sequence, plus, and qscore lines, and returns them as a list to be parsed'''
    # Note: Pete mentioned last time can just use readline in PS5
    header: str = file.readline().strip()
    sequence: str = file.readline().strip()
    plus: str = file.readline().strip()
    q_score: str = file.readline().strip()
    record: list = [header, sequence, plus, q_score]
    return record

def open_output_files(i_set):
    '''This function takes known index_set, opens output files created from the index-pair reads, and also opens files for hopped and unknown cases for R1 and R4'''
    fh_dict: dict = {}
    # setting fh indices for read 1 and read 4 as values 
    for index in index_seq:
        output_file_R1: str = index + "_R1.fq"
        output_file_R2: str = index + "_R4.fq"
        fh_matchedR1 = open(output_file_R1, "w")
        fh_matchedR2 = open(output_file_R2, "w")
        fh_dict[index] = [fh_matchedR1, fh_matchedR2]

    fh_hoppedR1 = open("hopped_R1.fq", "w")
    fh_hoppedR2 = open("hopped_R2.fq", "w")
    fh_dict["hopped"] = [fh_hoppedR1, fh_hoppedR2]

    fh_unknownR1 = open("unknown_R1.fq", "w")
    fh_unknownR2 = open("unknown_R2.fq", "w")
    fh_dict["unknown"] = [fh_unknownR1, fh_unknownR2]

    return fh_dict

output_file_dict: dict = open_output_files(index_seq)

with gzip.open(read1, "rt") as fh1, gzip.open(index1, "rt") as fh2, gzip.open(index2, "rt") as fh3, gzip.open(read2, "rt") as fh4:
    while True:
        r_counter += 1 #increment record counter 
        
        # get records in list form to parse through
        r1_list: list = create_record(fh1)
        i1_list: list = create_record(fh2)
        i2_list: list = create_record(fh3)
        r2_list: list = create_record(fh4)

        if r1_list[0] == "": #EOF 
            break
        
        rc_index2: str = rev_complement(i2_list[1])
        # print("Reverse complement", rc_index2)

        # print("R1\n", r1_list)
        # print("I1\n", i1_list)
        # print("I2\n", i2_list)

        r1_list[0], index_pair = append_to_header(r1_list[0], i1_list[1], rc_index2)
        r2_list[0], index_pair = append_to_header(r2_list[0], i1_list[1], rc_index2)

        
        # print("R2 changed header\n", r2_list, "\n")
                # if N in R2 (index1) or N in R3 record or R2 low quality or R3 low quality:
        if "N" in i1_list[1] or "N" in i2_list[1] or qual_check(i1_list[3])==False or qual_check(i2_list[3])==False:
            num_unknown += 1 # increment unknown_in 
            # write to file for read 1; use the dictionary to access file handle
            # Unknown R1 write header, sequence, plus, quality scores
            output_file_dict["unknown"][0].write(r1_list[0] + "\n" + r1_list[1] + "\n" + r1_list[2] + "\n" + r1_list[3] + "\n")
            # Unknown R4 write header, sequence, plus, quality scores
            output_file_dict["unknown"][1].write(r2_list[0] + "\n" + r2_list[1] + "\n" + r2_list[2] + "\n" + r1_list[3] + "\n")

        # elif index 1 or 2 not in index_seq (known indices)
        elif i1_list[1] not in index_seq or rc_index2 not in index_seq:
            num_unknown += 1 # increment unknown
            # write to file ...
            # Unknown R1 write header, sequence, plus, quality scores
            output_file_dict["unknown"][0].write(r1_list[0] + "\n" + r1_list[1] + "\n" + r1_list[2] + "\n" + r1_list[3] + "\n")
            # Unknown R4 write header, sequence, plus, quality scores
            output_file_dict["unknown"][1].write(r2_list[0] + "\n" + r2_list[1] + "\n" + r2_list[2] + "\n" + r1_list[3] + "\n")

        # elif R2 == R3 
        elif i1_list[1] == rc_index2:
            num_matched += 1 # increment matched counter
            # write to file
            if index_pair in matched_in:
                matched_in[index_pair] += 1
            else:
                matched_in[index_pair] = 1
            # Get the barcode/index sequence [1] for R1 [0] in dictionary, write header, sequence, plus, quality scores
            output_file_dict[i1_list[1]][0].write(r1_list[0] + "\n" + r1_list[1] + "\n" + r1_list[2] + "\n" + r1_list[3] + "\n")
            # Get the barcode/index sequence [1] for R4 [1] in dictionary, write header, sequence, plus, quality scores
            output_file_dict[i1_list[1]][1].write(r1_list[0] + "\n" + r1_list[1] + "\n" + r1_list[2] + "\n" + r1_list[3] + "\n")
        
        # elif R2 != R3
        elif i1_list != rc_index2:
            num_unmatched += 1 # unmatched, increment hopped
            # write to hopped
            if index_pair in hopped_in:
                hopped_in[index_pair] += 1
            else:
                hopped_in[index_pair] = 1
            # Hopped R1 write header, sequence, plus, quality scores
            output_file_dict["hopped"][0].write(r1_list[0] + "\n" + r1_list[1] + "\n" + r1_list[2] + "\n" + r1_list[3] + "\n")
            # Hopped R4 write header, sequence, plus, quality scores
            output_file_dict["hopped"][1].write(r2_list[0] + "\n" + r2_list[1] + "\n" + r2_list[2] + "\n" + r2_list[3] + "\n")

# Summary of Outputs to Terminal

total_counts: int = num_matched + num_unknown + num_unmatched
# Printing all counts
print(f'Total Records: {r_counter}')
print(f'Total Dual-Matched: {num_matched}')
print(f'Total Unknown: {num_unknown}')
print(f'Total Index Hopped: {num_unmatched}')
print(f'Total Counts: {total_counts}')

# Print Percentages of Matched, Unmatched/Hopped, Unknown
percent_matched: float = (num_matched/total_counts) * 100
print(f'Percent Matched: {percent_matched}')

percent_hopped: float = (num_unmatched/total_counts) * 100
print(f'Percent Hopped: {percent_hopped}')

percent_unknown: float = (num_unknown/total_counts) * 100
print(f'Percent Unknown: {percent_unknown}')

# Possible Index Pairs 
print("Possible Indexed Matched Pairs")
for i in matched_in:
    print(f'{i}\t{matched_in[i]}\tPercent:{(matched_in[i]/total_counts)*100}')
for j in hopped_in:
    print(f'{j}\t{hopped_in[j]}\tPercent:{(hopped_in[j]/total_counts)*100}')

# Summary of Outputs to File 
with open("Statistics_Summary.md", "wt") as stat_sum:
    stat_sum.write(f'---Summary Table---\n')
    stat_sum.write(f'Total Dual-matched: {num_matched}\nTotal Index Hopped: {num_unmatched}\nTotal Unknown: {num_unknown}\n')
    stat_sum.write(f'Total Index Hits {total_counts}\n')
    stat_sum.write(f'---Percentages---\n')
    stat_sum.write(f'Percent Matched: {percent_matched}\nPercent Hopped: {percent_hopped}\nPercent Unknown: {percent_unknown}\n')
    stat_sum.write(f'---Possible Indexed Matched Pairs:---\n')
    stat_sum.write("Dual-Matched Pairs")
    for k in matched_in:
        stat_sum.write(f'{k}\t{matched_in[k]}\tPercent:{(matched_in[k]/total_counts)*100}%')
    stat_sum.write(f'\nIndex-Hopped Pairs')
    for l in hopped_in:
        stat_sum.write(f'{l}\t{hopped_in[l]}\tPercent:{(hopped_in[l]/total_counts)*100}%')

# Plotting Matched Index-Pairs Distributions

x = matched_in.keys()
y = matched_in.values()
plt.xlabel("Index_pair")
plt.ylabel("Num Occurrences")
plt.title("Distribution of Matched Index Pairs from Demultiplexing")
plt.bar(x,y)
plt.savefig("Matched_Distribution.png")

#x2 = hopped_in.keys()
#y2 = hopped_in.values()



# Close 52 files (index matched, unknown, hopped/unmatched)
for key in output_file_dict:
    output_file_dict[key][0].close() # Closes R1 files
    output_file_dict[key][1].close() # Closes R4 (read 2) files

