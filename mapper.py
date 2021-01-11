##MAPPING ALG Project
"""
ALG Project : mapping without gap, SNPs detection
M2 IBI 2020-2021
Benjamin Blanc
Tiphaine Casy
"""
################IMPORT SECTION################
import pandas as pd
import argparse
import time
from collections import Counter

################PARSER SECTION################
parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="Path of genome file (.fasta)", type = str, default=argparse.SUPPRESS)
parser.add_argument("--index", help="Path of the file (.dp)",type = str, default=argparse.SUPPRESS)
parser.add_argument("--reads", help="Path of reads file (.fasta)", type = str, default=argparse.SUPPRESS)
parser.add_argument("--k", help="Size of windows for k-mers", type = int, default=argparse.SUPPRESS)
parser.add_argument("--max_hamming", help="Maximal number of substitution", type = int, default=argparse.SUPPRESS)
parser.add_argument("--min_abundance", help="Maximal abundance of variant", type = int, default=argparse.SUPPRESS)
parser.add_argument("--out", help="Path of out file (.vcf)", type = str, default=argparse.SUPPRESS)
args = parser.parse_args()

################TIME COUNT################
start = time.time()

################BWT INFORMATION AND RESEARCH################
def get_N(BWT):
    """
    Function used to count all A,T,C and G present in the reference sequence

    :param sequence: reference sequence with "$" at the end (i.e : "READ AND KEEPBACK SEQUENCE OF INPUT FILE" part)
    :type sequence: string

    :return:    For each nucleotide return its count
    :rtype:     dictionnary -> key : nucleotide (A,T,C or G) ; value : count
    """

    n = {"$": 0, "A": 0, "C":0, "G": 0, "T":0} # key = letter, value = number of occurrences
    for letter in BWT:
        n[letter] += 1
    return n

def R_table(BWT) :
    """
    Function used to indicate the rank of the BWT nucleotide (so, the last nucleotide of Burrows Wheeler sequence)
    example : if the first Burrows Wheeler sequence (so index "0") is "C1" this function returns "1"

    :param index: index of the nucleotide in the F list
    :param BWT: list of the nucleotide of BWT list

    :type index: int
    :type BWT: list of string

    :return:   list of dictionaries with count of each nucleotides already seen
    :rtype:    list
    """
    n = {} # key = letter, value = number of occurrences
    r = [] # for each i: rank of the i^th value in bwt
    for letter in BWT:
        if letter not in n:
            n[letter] = 0
        n[letter] += 1
        r.append(n[letter])
    return r

def get_querry(pattern: str, bwt: str, n: {}, r: [], sa: [int]) :
    """
    Function used to research a querry into the Burrows Wheeler

    :param BWT: list of last necleotides of Burrows Wheeler sequences
    :param Q: querry (searching sequence)
    :param N: dictionary created in "get_F" function
    :param sa: list of the index of the row in Burrows Wheeler sequences before sorting

    :type BWT: list of string
    :type Q: string
    :type N: dictionnary -> key : nucleotide (A,T,C or G) ; value : count
    :type sa: list of int

    :return:   position of founded querry into the reference sequence
    :rtype:    2 conditions : if there is one position -> int
                              if there are many positions -> list of int
    """
    start = 0
    stop = len(bwt)-1
    # read the pattern from right to left
    for pos_pattern in range (len(pattern)-1,-1,-1):
        current_char = pattern[pos_pattern]
        new_start = get_down(bwt, current_char, start, stop)
        if new_start == -1:
            return []
            break

        else :
            new_stop = get_up(bwt, current_char, start, stop)
            start = left_first(bwt[new_start], r[new_start], n)
            stop = left_first(bwt[new_stop], r[new_stop], n)

    res = []
    for i in range(start, stop+1):
        res.append(sa[i])
    return res

def get_down(bwt: str, alpha: chr, start: int, stop: int) -> int:
    """
    Detects the first occurrence of alpha in bwt for i in [start, stop].

    From start go down in the bwt as long as bwt[line] != alpha and line <= stop
      - if bwt[line] == alpha, returns the corresponding line
      - if line > stop: returns -1
    """
    line = start
    while line <= stop:
        if bwt[line] == alpha:
            return line
        line += 1
    return -1

def get_up(bwt: str, alpha: chr, start: int, stop: int) -> int:
    line = stop
    while line >= start:
        if bwt[line] == alpha:
            return line
        line -= 1
    return -1

def reverse_transcript(sequence):
    '''
    Return the reverse transcript of a given sequence

    :param sequence: A sequence of nucelotids
    :type sequence: string

    :return:   reversed transcript sequence
    :rtype:    string
    '''
    new_sequence = ''
    for x in reversed(range(len(sequence))):
        if "G" == sequence[x] :
            new_sequence += "C"
        if "A" == sequence[x]:
            new_sequence += "T"
        if "C" == sequence[x]:
            new_sequence+="G"
        if "T" == sequence[x]:
            new_sequence+="A"
    return(new_sequence)

def left_first(alpha: chr, k: int, n: {}) -> int:
    assert k<= n[alpha], f"Cannot ask for the {k}^th {alpha}, it does not exist"
    if alpha == "$": return 0
    if alpha == "A": return n["$"] + k - 1
    if alpha == "C": return n["$"] + n["A"] + k - 1
    if alpha == "G": return n["$"] + n["A"] + n["C"] + k - 1
    if alpha == "T": return n["$"] + n["A"] + n["C"] + n["G"] + k - 1
    raise ValueError(f"Character {alpha} not in the bwt")

################OPEN READS FILE################
def search_querry(reads, k_mer, index, N, R, max_hamming) :
    """
        This function searches correspondance(s) of kmers in a sequence thanks to its index.

        :param reads: A fasta file containing sequences of nucleotids
        :param k_mer: The length used for the search of kmers
        :param index: The Burrows Wheeler index of a sequence, in format .dp
        :type reads: fasta file
        :type kmer: int
        :type index: dp file
        :return: Return, for every kmer in given reads, the position in a sequence where it aligns, the number of the kmer in the read and the read direction.
                 Return the list of the reads for each sense as well.
        :rtype: dictionary

        :Example:

        An example of output is :

        {('TCTGA', 0): [[301, 500], 0, '+'], ('CTGAT', 1): [501, 1, '+']}

        The key of the dictionary includes the kmer and the position of this kmer in the read.

        The values for each key respectively includes :
        - one or several positions in the genom correspondingto the kmer
        - the position of this kmer in the read
        - the sens of the strand where the correspondance is found, + for sense and - for antisense

    """
    print("BEGIN SEARCH")
    list_querry = []
    read_information = {}
    for read_lines in reads : #Read the file in both senses
        read_lines = str(read_lines).replace('\n','')
        if '>' not in read_lines : #Indicates which read is on which strand.
            reverse_read_lines = reverse_transcript(str(read_lines))
            kmer_sai = {}
            k_mer_modified = k_mer
            if "+" not in read_information.keys() and "-" not in read_information.keys():
                read_information['+'] = [str(read_lines)]
                read_information['-'] = [str(reverse_read_lines)]
            else :
                read_information['+'].append(str(read_lines))
                read_information['-'].append(str(reverse_read_lines))

            for x in range(0, len(read_lines)-k_mer_modified+1): #Search correspondance in the genome
                sai_querry = get_querry(str(read_lines[x:k_mer_modified]),index['BWT'], N, R, index['SA[i]'])
                if sai_querry : #Verify if the list is full Find a new correspondance on the sense strand and add it in the dictionary
                    kmer_sai[str(read_lines[x:k_mer_modified]),x] = [sai_querry,x, "+"] ##(sai_querry,x)

                else : #Find a new correspondance on the antisense strand and add it in the dictionary
                    sai_querry = get_querry(str(reverse_read_lines[x:k_mer_modified]),index['BWT'],N, R, index['SA[i]'])
                    kmer_sai[str(reverse_read_lines[x:k_mer_modified]),x] = [sai_querry,x, "-"] ##(sai_querry,x)

                k_mer_modified +=1
            list_querry.append(kmer_sai)
    print("END OF SEARCH")

    return(list_querry,read_information)



def comparison(or_sequence, reads, start_comparison, max_hamming) :
    '''
        Align two sequences from a given position and return the number of substition and the starting alignment position

        :param or_sequence: A sequence of nucleotide
        :param reads: A fasta file containing sequences of nucleotids
        :param start_comparison: The position of the beginning of alignment
        :param max_substitution: The maximum of differences allowed between both sequences
        :type or_sequence: string
        :type reads: fasta file
        :type start_comparison: int
        :type max_substitution: int
        :return: The number of substitutions and the position of the beginning of alignement
        :rtype: list
    '''

    index_subs = []

    comparison_changed = start_comparison
    if start_comparison <= len(or_sequence)-len(reads) :
        substitution = 0
        for x in range(len(reads)): #For each read, save the substitutions
            if or_sequence[comparison_changed] is not reads[x] and substitution <= max_hamming:
                substitution +=1
                result = (or_sequence[comparison_changed], comparison_changed, reads[x])
                index_subs.append(result)
            comparison_changed+=1

    if len(index_subs) <= max_hamming:
        p = [start_comparison,index_subs]
    else :
        p = ['','']
    return(p)

def seed_and_extend(sequence, querry_found, max_hamming, min_abundance, output) :
    """
        Displays the read ID and its best mapping position and the number of substitions observed

        :param sequence: A sequence of nucleotide
        :param querry_found: The output of the function search_querry()
        :param max_hamming: The maximum number of substitutions allowed
        :type sequence: string
        :type querry_found: dictionary
        :type max_hamming: int
    """
    print("BEGIN SEED AND EXTEND")
    sequence = sequence[:-1]
    sai_of_kmer = querry_found[0]
    reads = querry_found[-1]
    substituion_info = []
    #print(sai_of_kmer)
    for x in range(0,len(sai_of_kmer)) : ##For each read
        comparison_result = {}
        list_of_position = []
        for y in range(0,len(sai_of_kmer[x])): ##For each kmer of a read
            sai_values = list(sai_of_kmer[x].values())[y][0]

            for i in sai_values :
                if i != '':
                    index_kmer = list(sai_of_kmer[x].values())[y][1]
                    k_mer_sens = list(sai_of_kmer[x].values())[y][-1]
                    start_comparison = int(i-index_kmer)
                    read_good_sense = reads[k_mer_sens][x]
                    if start_comparison >= 0 and start_comparison not in list_of_position:
                        list_of_position.append(start_comparison)
                    elif start_comparison in list_of_position  : break

            for pos in list_of_position :
                results = comparison(sequence, read_good_sense, pos, max_hamming)
                if (results[0]) not in comparison_result.keys() and results[0] != '': #If the number of substitutions is not already in the dictionary, create a new key
                    comparison_result[(results[0])] = [(results[1],k_mer_sens)]

        for key in comparison_result.keys():
            for i in comparison_result[key][0][0] :
                substituion_info.append(i)

    vcf_creation = Counter(substituion_info)
    position_subsitution = []
    original_nucleotide = []
    reads_nucleotide = []
    number_substitution = []

    for key in vcf_creation.keys():
        if int(vcf_creation[key]) >= int(min_abundance) :
            position_subsitution.append(key[1])
            original_nucleotide.append(key[0])
            reads_nucleotide.append(key[-1])
            number_substitution.append(vcf_creation[key])
        else :
            break
    d = {'POS' : position_subsitution, 'REF' : original_nucleotide, 'ALT' : reads_nucleotide, 'ABUNDANCE' : number_substitution}
    df = pd.DataFrame(data = d)
    df = df.sort_values('POS')
    df.to_csv(str(output), index = False, encoding= 'utf-8', mode = 'a', header = False, sep='\t')
    print("END")

################CHECK IF PARSER IS FULL################


if format(args) != 'Namespace()':
    ref = open(str(args.ref), 'r')
    index = pd.read_csv(str(args.index))
    reads = open(str(args.reads), 'r')
    k_mer = args.k
    max_hamming = args.max_hamming
    min_abundance = args.min_abundance
    output = args.out

    information_file = open(str(output), 'w')
    information_file.write('##REFERENCE : ' + str('reference.fasta') + '\n')
    information_file.write('##READS : ' + str('reads.fasta') + '\n')
    information_file.write('##K : ' + str(k_mer) + '\n')
    information_file.write('##MAX_SUBS : ' + str(max_hamming) + '\n')
    information_file.write('##MIN_ABUNDANCE : ' + str(min_abundance) + '\n')
    information_file.write('##SUBSTITUTION INFORMATION \n')
    information_file.write('#POS' + '\t' + str('REF') + '\t' + str('ALT') + '\t' + str('ABUNDANCE INFORMATION') + '\n')
    information_file.close()

    ################READ AND KEEPBACK SEQUENCE OF INPUT FILE################
    sequence = ''
    for line in ref :
        line = str(line).replace('\n','')
        if '>' not in line :
            sequence = str(line) + '$'
    seed_and_extend(sequence, search_querry(reads, k_mer, index,get_N(index["BWT"]),R_table(index["BWT"]), max_hamming), max_hamming, min_abundance, output)
else :
    print("Obligation to inform all argues to search substitutions \n")
    print("For some help write : \'mapper.py --help\' in control terminal")

################TIME COUNT################
end = time.time()

time_to_analyse = end-start
print("TIME FOR ANALYSE : " + str(round(time_to_analyse,3)) + " secondes")
