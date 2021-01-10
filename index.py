"""
ALG Project : mapping without gap, SNPs detection
M2 IBI 2020-2021
Benjamin Blanc
Tiphaine Casy
"""
################IMPORT SECTION################
import tools_karkkainen_sanders as tks
import pandas as pd
import argparse
import time

################PARSER SECTION################
parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="Path of genome file (.fasta)", type = str, default=argparse.SUPPRESS)
parser.add_argument("--out", help="Path of the results file (.dp)", type = str, default=argparse.SUPPRESS)
args = parser.parse_args()


################TIME COUNT################
start = time.time()

################F,i,BWT KEEP BACK################
def get_BWT(sequence,sa):
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
    bwt = []
    for i in range(len(sa)):
        bwt.append(sequence[sa[i]-1])
    return bwt
################CHECK IF PARSER IS FULL################
if format(args) != 'Namespace()':
    ref = open(str(args.ref), 'r')

    ################READ AND KEEPBACK SEQUENCE OF INPUT FILE################
    """
    To use the Burrows Wheeler algorithm we need to add a "$" at the end of sequence to be usable by the other functions
    """
    sequence = ''
    for line in ref :
        line = str(line).replace('\n','') ##add to delete line break
        if '>' not in line :
            sequence = str(line)+"$"


    ################SA[i] KEEP BACK################
    sa = tks.simple_kark_sort(sequence) ##keep back of SA[i] calculated thanks to tools_karkkainen_sanders
    get_BWT(sequence,sa)

    ################CREATION OF OUT FILE################
    """
    Save the index in a file as dataframe format, path indication in argues at the beginning
    """

    d = {'SA[i]' : sa,'BWT' : get_BWT(sequence,sa)}
    df = pd.DataFrame(data = d)
    df.to_csv(str(args.out), encoding= 'utf-8', index=False, mode = 'w', header = True)

else :
    print("Obligation to inform all argues to search substitutions \n")
    print("For some help write : \'mapper.py --help\' in control terminal")




################TIME COUNT################
end = time.time()

time_to_analyse = end-start
print("TIME FOR ANALYSE : " + str(round(time_to_analyse,3)) + " secondes")
