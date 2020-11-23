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
parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="Path of genome file (fasta)")
parser.add_argument("--out", help="Path of the results file (.dp)")
args = parser.parse_args()

################READ AND KEEPBACK SEQUENCE OF INPUT FILE################
file = open(str(args.ref), 'r')
data = file.readlines()
file.close()
sequence = []
for line in data :
    line=line.replace('\n','')
    if '>' not in line :
        sequence.append(line)
sequence = ''.join(sequence[:-1]) + "$"  ##fasta file have a particular format

################SA[i] KEEP BACK################
sa = tks.simple_kark_sort(sequence) ##keep back of SA[i] calculated thanks to tools_karkkainen_sanders

################F,i,BWT KEEP BACK################
def get_BWT(sequence):
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
    return BWT, F, BWT_list
results_BWT = get_BWT(sequence)

################CREATION OF OUT FILE################
### writing into a panda frame all information for reference sequence
### and save it into a '.dp' file 
d = {'SA[i]' : sa, 'F' : results_BWT[1], 'Sequence_BWT' : results_BWT[-1], 'BWT' : results_BWT[0]}
df = pd.DataFrame(data = d)
df.to_csv(str(args.out), encoding= 'utf-8', index=False, mode = 'w', header = True)
