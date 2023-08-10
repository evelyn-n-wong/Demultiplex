# Author: Evelyn Wong evew@uoregon.edu, wonge29@gmail.com
#!/usr/bin/env python 
# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

#Last update: 2023-August-9 __version__0.5

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33 

def qual_score(phred_score: str) -> float:
    '''This function takes in a string, parses through the characters
      in the string, adds the the individual Q-score, and returns the average Q-score'''
    length = len(phred_score)
    sum: int = 0
    for char in phred_score: #goes through individual characters in the string
        sum = sum + convert_phred(char)
    return sum/length

def validate_base_seq(seq: str, RNAflag: bool = False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    #set_seq = set(seq)
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def calc_median(sorted_sl: list) -> list:
    '''This function takes in a one dimensional list and returns one median value'''
    median: int
    s_length = len(sorted_sl) // 2
    if len(sorted_sl) % 2 == 0:
        median = ((sorted_sl[s_length-1] + sorted_sl[s_length])/2)
    else:
        median = (sorted_sl[s_length])
    return median

def gc_content(sequence):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(sequence) == True, "Valid Base Sequence inputed"
    sum_GC = sequence.count("G") + sequence.count("C")
    length_DNA = len(sequence)
    return (sum_GC/length_DNA)

def oneline_fasta(input_file):
    '''This function takes in a multiline FASTA file and returns a one-line fasta'''
    with open(input_file, "r") as fh, open({input_file}_oneline.fa, "w") as f_out:
#    fh.readline() # gets the first line. Also, everytime this is called, subsequent lines are read sequentially.
        x: int = 0 #make sure it's the first line else
        seq: str = ""
        for line in fh:
            if not line.startswith(">"):
                seq += line.strip('\n')
            if line.startswith(">"):
                if x == 0:
                    f_out.write(line) #header
                else:
                    f_out.write(seq + "\n")
                    f_out.write(line)
                    seq = ""
            x+=1
        f_out.write(seq) # gets last line  

print("working in", __name__, "namespace") #when 

if __name__ == "__main__":
    #tests for validate_base_seq()
    assert validate_base_seq("AAA") == True
    assert validate_base_seq("xxx") == False
    print("Sequences passed Validation of Bases tests")
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("GC content correctly calculated!")
