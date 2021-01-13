##MAPPING ALG Project
"""
ALG Project : mapping without gap, SNPs detection
M2 IBI 2020-2021
Benjamin Blanc ~ Tiphaine Casy
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
def get_N(bwt):
    """
    Function used to count all A,T,C and G present in the reference sequence

    :param bwt: list of Burrows Wheeler Transform of origin sequence
    :type bwt: list of string

    :return:    For each nucleotide return its count
    :rtype:     dictionnary -> key : nucleotide (A,T,C or G) ; value : count
    """

    n = {"$": 0, "A": 0, "C":0, "G": 0, "T":0} # key = letter, value = number of occurrences
    for letter in BWT:
        n[letter] += 1
    return n

def R_table(bwt) :
    """
    Function used to indicate the rank of the BWT nucleotide
    example : if the first Burrows Wheeler sequence (so index "0") is "C1" this function returns "1"

    :param bwt: list of Burrows Wheeler Transform of origin sequence
    :type bwt: list of string

    :return:   list of dictionary with count of each nucleotides already seen
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

def get_querry(pattern, bwt, n, r, sa) :
    """
    Function used to research a querry into the Burrows Wheeler

    :param pattern: querry research
    :type pattern: str

    :param bwt: list of Burrows Wheeler Transform of origin sequence
    :type bwt: list of string

    :param n: dictionnary of count of each type of nucleotides in the sequence
    :type n: dictionnary -> key : nucleotide (A,T,C or G) ; value : count

    :param r: list of dictinonary with count of each nucleotides already seen
    :type r: list

    :param sa: list of the index of the row in Burrows Wheeler sequences before sorting
    :type sa: list of int

    :return: list of position of founded querry into the reference sequence, if querry isn't found the list is empty
    :rtype: list of int
    """
    start = 0
    stop = len(bwt)-1
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

def get_down(bwt, alpha, start, stop):
    """
    Detects the first occurrence of alpha in bwt for i in [start, stop].

    From start go down in the bwt as long as bwt[line] != alpha and line <= stop
      - if bwt[line] == alpha, returns the corresponding line
      - if line > stop: returns -1

    :param bwt: list of Burrows Wheeler Transform of origin sequence
    :type bwt: list of string

    :param alpha: nucleotide wanted
    :type alpha: str

    :param start: position of first occurence of alpha in F column
    :type start: int

    :param stop: position of last occurence of alpha in F column
    :type stop: int

    :return: position of alpha wanted
    :rtype: int
    """
    line = start
    while line <= stop:
        if bwt[line] == alpha:
            return line
        line += 1
    return -1

def get_up(bwt, alpha, start, stop):
    """
    Detects the last occurrence of alpha in bwt for i in [start, stop].

    From start go down in the bwt as long as bwt[line] != alpha and line <= stop
      - if bwt[line] == alpha, returns the corresponding line
      - if line > stop: returns -1

    :param bwt: list of Burrows Wheeler Transform of origin sequence
    :type bwt: list of string

    :param alpha: nucleotide wanted
    :type alpha: str

    :param start: position of first occurence of alpha in F column
    :type start: int

    :param stop: position of last occurence of alpha in F column
    :type stop: int

    :return: position of alpha wanted
    :rtype: int
    """
    line = stop
    while line >= start:
        if bwt[line] == alpha:
            return line
        line -= 1
    return -1

def reverse_transcript(sequence):
    '''
    Return the reverse transcript of a given sequence

    :param sequence: sequence of nucleotide (here : sequence of reads)
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

def left_first(alpha, k, n):
    """
    Return the number of the row of first value of lexicographic ordered sequences for alpha nucleotide

    :param alpha: nucleotide wanted
    :type alpha: str

    :param k: rank of the letter alpha in the bwt
    :type k: int

    :param n: dictionnary of count of each type of nucleotides in the sequence
    :type n: dictionnary -> key : nucleotide (A,T,C or G) ; value : count

    :return: number of the row of alpha
    :rtype: int
    """
    assert k<= n[alpha], f"Cannot ask for the {k}^th {alpha}, it does not exist"
    if alpha == "$": return 0
    if alpha == "A": return n["$"] + k - 1
    if alpha == "C": return n["$"] + n["A"] + k - 1
    if alpha == "G": return n["$"] + n["A"] + n["C"] + k - 1
    if alpha == "T": return n["$"] + n["A"] + n["C"] + n["G"] + k - 1
    raise ValueError(f"Character {alpha} not in the bwt")

################OPEN READS FILE################
def search_querry(reads, k_mer, index, N, R) :
    """
    This function searches correspondance(s) of kmers in a sequence thanks to its index.

    :param reads: A fasta file containing sequences of nucleotids
    :type reads: fasta file

    :param k_mer: The length used for the search of kmers
    :type kmer: int

    :param index: The Burrows Wheeler index of a sequence, in format .dp
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
            if "+" not in read_information.keys() and "-" not in read_information.keys():
                read_information['+'] = [str(read_lines)]
                read_information['-'] = [str(reverse_read_lines)]
            else :
                read_information['+'].append(str(read_lines))
                read_information['-'].append(str(reverse_read_lines))

            kmer_list = k_mer_creation(str(read_lines), int(k_mer))
            kmer_reversed_list = k_mer_creation(str(reverse_read_lines), int(k_mer))
            for x in range(0,len(kmer_list)) :
                sai_querry = get_querry(kmer_list[x],index['BWT'], N, R, index['SA[i]'])
                if sai_querry : #Verify if the list is full Find a new correspondance on the sense strand and add it in the dictionary
                    kmer_sai[kmer_list[x],x] = [sai_querry,x, "+"] ##(sai_querry,x)
                else : #Find a new correspondance on the antisense strand and add it in the dictionary
                    sai_querry = get_querry(kmer_reversed_list[x],index['BWT'],N, R, index['SA[i]'])
                    kmer_sai[kmer_reversed_list[x],x] = [sai_querry,x, "-"] ##(sai_querry,x)

            list_querry.append(kmer_sai)
    print("END OF SEARCH")
    return(list_querry,read_information)

def k_mer_creation(sequence, kmer_size):
    """
    Create the list of k_mer of sequence depending of the kmer_size choosen

    :param sequence: sequence of nucleotide (here : sequence of reads)
    :type sequence: string

    :param kmer_size: Size of length of cut origin sequence
    :type kmer_size: int

    :return: list of cut sequence with length of kmer_size
    :rtype: list of str
    """
    k_mer = []
    k_mer_modified = kmer_size
    for x in range(0, len(sequence)-k_mer_modified+1):
        k_mer.append(str(sequence[x:k_mer_modified]))
        k_mer_modified += 1
    return(k_mer)

def comparison(or_sequence, reads, start_comparison, max_hamming) :
    '''
    Align two sequences from a given position and return the number of substition and the starting alignment position
    :param or_sequence: A sequence of nucleotide
    :type or_sequence: string

    :param reads: A fasta file containing sequences of nucleotids
    :type reads: fasta file

    :param start_comparison: The position of the beginning of alignment
    :type start_comparison: int

    :param max_substitution: The maximum of differences allowed between both sequences
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
    Search the best mapping position of reads on reference sequence and return all differences betwenn reads and reference sequences (SNPs)

    :param sequence: reference sequence
    :type sequence: reference sequence

    :param querry_found: The output of the function search_querry()
    :type querry_found: dictionary

    :param max_hamming: the maximum number of substition allowed
    :type max_hamming: int

    :param min_abundance: the minimum of repetition of SNP in the reference sequence
    :type min_abundance: int

    :param output: path of output file (.vcf) to write all information of mapping
    :type output: str
    """
    print("BEGIN SEED AND EXTEND")
    sequence = sequence[:-1]
    sai_of_kmer = querry_found[0]
    reads = querry_found[-1]
    substitution_info = []
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
                substitution_info.append(i)

    vcf_creation = Counter(substitution_info)
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

    ################READ AND KEEPBACK SEQUENCE OF INPUT FILE###############
    sequence = ''
    for line in ref :
        line = str(line).replace('\n','')
        if '>' not in line :
            sequence = str(line) + '$'

    information_file = open(str(output), 'w')
    information_file.write('##REFERENCE : ' + str(args.ref) + '\n')
    information_file.write('##READS : ' + str(args.reads) + '\n')
    information_file.write('##K : ' + str(k_mer) + '\n')
    information_file.write('##MAX_SUBS : ' + str(max_hamming) + '\n')
    information_file.write('##MIN_ABUNDANCE : ' + str(min_abundance) + '\n')
    information_file.write('##SUBSTITUTION INFORMATION \n')
    information_file.write('#POS' + '\t' + str('REF') + '\t' + str('ALT') + '\t' + str('ABUNDANCE INFORMATION') + '\n')
    information_file.close()

    seed_and_extend(sequence, search_querry(reads, k_mer, index,get_N(index["BWT"]),R_table(index["BWT"])), max_hamming, min_abundance, output)

else :
    print("Obligation to inform all argues to search substitutions \n")
    print("For some help write : \'mapper.py --help\' in control terminal")

################TIME COUNT################
end = time.time()

time_to_analyse = end-start
print("TIME FOR ANALYSE : " + str(round(time_to_analyse,3)) + " secondes")
