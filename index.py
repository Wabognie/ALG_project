"""
ALG Project : mapping without gap, SNPs detection
M2 IBI 2020-2021
Benjamin Blanc
Tiphaine Casy
"""
################IMPORT SECTION################
import tools_karkkainen_sanders as tks
import pandas as pd
import sys
import argparse

################PARSER SECTION################
"""
TO DO : try execpt pour obliger a mettre les arguments
"""

"""
    Try - except parser : force to enquire the 2 argues to run correctly the code
    "help" part : if user need some information to run code write "python3 index.py -help" so all informations are write on terminal
    To use argues and return the user input : "args.nameofargue" (ex : "args.ref" or "args.out" in this case)
"""
parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="Path of genome file (.fasta)")
parser.add_argument("--out", help="Path of the results file (.dp)")
args = parser.parse_args()

################READ AND KEEPBACK SEQUENCE OF INPUT FILE################
"""
    To use the Burrows Wheeler algorithm we need to add a "$" at the end of sequence to be usable by the other functions
"""
ref = open(str(args.ref), 'r')
sequence = ''
for line in ref :
    line = str(line).replace('\n','') ##add to delete line break
    if '>' not in line :
        sequence = str(line) + '$'

################SA[i] KEEP BACK################
sa = tks.simple_kark_sort(sequence) ##keep back of SA[i] calculated thanks to tools_karkkainen_sanders

################F,i,BWT KEEP BACK################
def get_BWT(sequence):
    """
    Function used to keep back sorted Burrows Wheeler sequences

    :param sequence: reference sequence with "$" at the end (i.e : "READ AND KEEPBACK SEQUENCE OF INPUT FILE" part)
    :type sequence: string

    :return:    'BWT' corresponding to the last caractere of sorted Burrows Wheeler sequences
                'F' corresponding to the first caractere of sorted Burrows Wheeler sequences
    :rtype:     list
                list
    """
    ### function used to keep back sorted BW sequences into a list('BWT_list')
    ### first letter of BW sequences ('F') and last one ('BWT')
    BWT_list = []
    BWT_list.append(sequence)
    for i in range(len(sequence)-1) :
        sequence = sequence[1:]+sequence[0]
        if str(sequence) not in BWT_list :
            BWT_list.append(sequence)

    BWT_list = sorted(BWT_list)
    BWT = []
    F = []
    for i in BWT_list :
        BWT.append(i[-1])
        F.append(i[0])
    return BWT, F
results_BWT = get_BWT(sequence)

################CREATION OF OUT FILE################
"""
Save the index in a file as dataframe format, path indication in argues at the beginning
"""

d = {'SA[i]' : sa, 'F' : results_BWT[1], 'BWT' : results_BWT[0]}
df = pd.DataFrame(data = d)
df.to_csv(str(args.out), encoding= 'utf-8', index=False, mode = 'w', header = True)
