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
def get_N(sequence):
    """
    Function used to count all A,T,C and G present in the reference sequence

    :param sequence: reference sequence with "$" at the end (i.e : "READ AND KEEPBACK SEQUENCE OF INPUT FILE" part)
    :type sequence: string

    :return:    For each nucleotide return its count
    :rtype:     dictionnary -> key : nucleotide (A,T,C or G) ; value : count
    """
    N = {}
    for i in sequence :
        if i not in N.keys():
            N[i]=1
        else :
            N[i]+=1
    return N

def LF(alpha, k,N):
    """
    Function used to attribute a rank for each ordered nucleotides (first value of Burrows Wheeler sequence, 'F'), starting with "$"

    :param alpha: reading nucleotide
    :param k: number of line in Burrows Wheeler index
    :param N: dictionary created in "get_F" function

    :type alpha: string
    :type k: int
    :type N: dictionnary -> key : nucleotide (A,T,C or G) ; value : count


    :return:    number of each sorted first nucleotides of Burrows Wheeler sequcence ('F')
    :rtype:     int
    """

    index = 0
    if alpha == "$":
        return(index)

    if alpha == "A" :
        index = int(N["$"])+k-1
        return(index)

    if alpha == "C" :
        index = int(N["$"])+int(N["A"])+k-1
        return(index)

    if alpha == "G" :
        index = int(N["$"])+int(N["A"])+int(N["C"])+k-1
        return(index)

    if alpha == "T" :
        index = int(N["$"])+int(N["A"])+int(N["C"])+int(N["G"])+k-1
        return(index)

def R_table(index,BWT) :
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
    R = []
    N = {}
    for letter in BWT :
        if letter not in N :
            N[letter] = 0
        N[letter] += 1
        R.append(N[letter])

    return R[index]

def get_querry(BWT, Q, N, sa) :
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

    last_char = Q[-1]
    i = LF(last_char,1,N)
    j = LF(last_char, N[last_char], N)
    i_min = -1
    j_max = -1

    breaker = False
    for x in reversed(range(len(Q)-1)):
        index = []
        for l in range(i, j+1):
            if str(BWT[l]) == str(Q[x]) :
                index.append(l)
                i_min = min(index)
                j_max = max(index)

        if len(index) == 0:
            breaker = True
        if i_min >=0 and j_max >=0:
            i = LF(BWT[i_min],R_table(i_min, BWT),N)
            j = LF(BWT[j_max],R_table(j_max, BWT),N)

        else :
            breaker = True

    found_sai = []
    if breaker :
        return('')
    if int(i) == int(j) :
        return(sa[i])
    for t in range(i,j+1) :
        found_sai.append(sa[t])

    return(sorted(found_sai))

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

################OPEN READS FILE################
def search_querry(reads, k_mer, index) :
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
    end = time.time()
    time_to_analyse = end-start
    print("TIME FOR ANALYSE BEGIN SEARCH: " + str(round(time_to_analyse,3)) + " secondes")
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
                sai_querry = get_querry(index['BWT'], str(read_lines[x:k_mer_modified]), get_N(sequence), index['SA[i]'])

                if sai_querry != '': #Find a new correspondance on the sense strand and add it in the dictionary
                    kmer_sai[str(read_lines[x:k_mer_modified]),x] = [sai_querry,x, "+"] ##(sai_querry,x)

                else : #Find a new correspondance on the antisense strand and add it in the dictionary
                    sai_querry = get_querry(index['BWT'],str(reverse_read_lines[x:k_mer_modified]), get_N(sequence), index['SA[i]'])
                    kmer_sai[str(reverse_read_lines[x:k_mer_modified]),x] = [sai_querry,x, "-"] ##(sai_querry,x)

                k_mer_modified +=1
            list_querry.append(kmer_sai)

    end = time.time()
    time_to_analyse = end-start
    print("TIME FOR ANALYSE END SEARCH: " + str(round(time_to_analyse,3)) + " secondes")    

    return(list_querry,read_information)



def comparison(or_sequence, reads, start_comparison, max_substitution, index_k_mer) :
    '''
        Align two sequences from a given position and return the number of substition and the starting alignment position

        :param or_sequence: A sequence of nucleotide
        :param reads: A fasta file containing sequences of nucleotids
        :param start_comparison: The position of the beginning of alignment
        :param max_substitution: The maximum of differences allowed between both sequences
        :param index_k_mer: The position of a kmer, in the read
        :type or_sequence: string
        :type reads: fasta file
        :type start_comparison: int
        :type max_substitution: int
        :type index_k_mer: int
        :return: The number of substitutions and the position of the beginning of alignement
        :rtype: list
    '''


    index_subs = []

    comparison_changed = start_comparison
    if start_comparison <= len(or_sequence)-len(reads) :
        substitution = 0
        for x in range(len(reads)): #For each read, save the substitutions
            if or_sequence[comparison_changed] != reads[x]:
                substitution +=1
                result = (or_sequence[comparison_changed], comparison_changed, reads[x])
                index_subs.append(result)
            comparison_changed+=1

        if len(index_subs) <= max_substitution:
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
    end = time.time()

    time_to_analyse = end-start
    print("TIME FOR ANALYSE BEGIN SEED: " + str(round(time_to_analyse,3)) + " secondes")
    sequence = sequence[:-1]
    sai_of_kmer = querry_found[0]
    reads = querry_found[-1]

    substituion_info = []

    for x in range(0,len(sai_of_kmer)) : ##For each read
        comparison_result = {}
        list_of_position = []
        for y in range(0,len(sai_of_kmer[x])): ##For each kmer of a read
            sai_values = list(sai_of_kmer[x].values())[y][0]
            index_kmer = list(sai_of_kmer[x].values())[y][1]
            k_mer_sens = list(sai_of_kmer[x].values())[y][-1]

            read_good_sense = reads[k_mer_sens][x]
            if isinstance(sai_values, list): ##verify if there is several positions of alignment
                for i in sai_values : #For each correspondance
                    if i != '':
                        start_comparison = int(i-index_kmer)
                        if start_comparison >= 0 and start_comparison not in list_of_position:
                            list_of_position.append(start_comparison)
            else : # One correspondance only
                if sai_values != '' :
                    start_comparison = int(sai_values-index_kmer)
                    if start_comparison >= 0 and start_comparison not in list_of_position:
                        list_of_position.append(start_comparison)

            for pos in list_of_position :
                results = comparison(sequence, read_good_sense, pos, max_hamming, index_kmer)
                if (results[0]) not in comparison_result.keys(): #If the number of substitutions is not already in the dictionary, create a new key
                    comparison_result[(results[0])] = [(results[1],k_mer_sens)]
                else : #Else, add the results to the existing key
                    comparison_result[(results[0])].append((results[1],k_mer_sens))

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
    d = {'POS' : position_subsitution, 'REF' : original_nucleotide, 'ALT' : reads_nucleotide, 'ABUNDANCE' : number_substitution}
    df = pd.DataFrame(data = d)
    df = df.sort_values('POS')
    df.to_csv(str(output), index = False, encoding= 'utf-8', mode = 'a', header = True, sep='\t')


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
    information_file.write('REFERENCE : ' + str('reference.fasta') + '\n')
    information_file.write('READS : ' + str('reads.fasta') + '\n')
    information_file.write('K : ' + str(k_mer) + '\n')
    information_file.write('MAX_SUBS : ' + str(max_hamming) + '\n')
    information_file.write('MIN_ABUNDANCE : ' + str(min_abundance) + '\n')
    information_file.write('\n')
    information_file.write('SUBSTITUTION INFORMATION \n')
    information_file.close()

    ################READ AND KEEPBACK SEQUENCE OF INPUT FILE################
    sequence = ''
    for line in ref :
        line = str(line).replace('\n','')
        if '>' not in line :
            sequence = str(line) + '$'

    seed_and_extend(sequence, search_querry(reads, k_mer, index), max_hamming, min_abundance, output)
else :
    print("Obligation to inform all argues to search substitutions \n")
    print("For some help write : \'mapper.py --help\' in control terminal")

################TIME COUNT################
end = time.time()

time_to_analyse = end-start
print("TIME FOR ANALYSE : " + str(round(time_to_analyse,3)) + " secondes")
