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

################PARSER SECTION################
parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="Path of genome file (.fasta)")
parser.add_argument("--index", help="Path of the file (.dp)")
parser.add_argument("--reads", help="Path of reads file (.fasta)")
parser.add_argument("--k", help="Size of windows for k-mers")
parser.add_argument("--max_hamming", help="Maximal number of substitution")
parser.add_argument("--min_abundance", help="Maximal abundance of variant")
parser.add_argument("--out", help="Path of out file (.vcf)")
args = parser.parse_args()

################TIME COUNT################
start = time.time()

################OPEN FILES################
ref = open('./reference.fasta','r') ##args.ref
index = pd.read_csv('./index.dp') ##args.index
k_mer = 20 ##args.k
max_hamming = 4 ##args.max_hamming

################CREATE SEQUENCE WITH $################
sequence = ''
for line in ref :
    line = str(line).replace('\n','')
    if '>' not in line :
        sequence = str(line) + '$'

################BWT INFORMATION AND RESEARCH################
def get_N(sequence):
    N = {}
    for i in sequence :
        if i not in N.keys():
            N[i]=1
        else :
            N[i]+=1
    return N

def LF(alpha, k,N):
    index = 0
    if alpha == "$":
        ##print(index)
        return(index)

    if alpha == "A" :
        index = int(N["$"])+k-1
        #print(index)
        return(index)

    if alpha == "C" :
        index = int(N["$"])+int(N["A"])+k-1
        #print(index)
        return(index)

    if alpha == "G" :
        index = int(N["$"])+int(N["A"])+int(N["C"])+k-1
        #print(index)
        return(index)

    if alpha == "T" :
        index = int(N["$"])+int(N["A"])+int(N["C"])+int(N["G"])+k-1
        #print(index)
        return(index)

def R_table(index,BWT) :
    ##retourne l'index de la lettre dans la liste F : si 0 de BWT = C1 -> return de 1
    R = [] # Liste des rangs dans la BWT
    N = {} # Dictionnaire pour retenir les occurrences des caractères qu'on a déjà rencontré
    for letter in BWT :
        if letter not in N :
            N[letter] = 0
        N[letter] += 1
        R.append(N[letter])

    return R[index]

def get_querry(BWT, Q, N, sa) :
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
#get_querry(get_BWT(sequence)[0], "GG", get_N(sequence), sa)

def reverse_transcript(sequence):
    '''
        Return the reverse transcript of a given sequence
        :param sequence: A sequence of nucelotids
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
    list_querry = []
    read_information = {}
    for read_lines in reads : #Read the file in both senses
        read_lines = str(read_lines).replace('\n','')
        reverse_read_lines = reverse_transcript(str(read_lines))

        if '>' not in read_lines : #Indicates which read is on which strand.
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
                    #if sai_querry == '':
                        #kmer_sai[str(reverse_read_lines[x:k_mer_modified]),x] = ["NONE",x, "None"]

                k_mer_modified +=1
            list_querry.append(kmer_sai)

    return(list_querry,read_information)
querry_found = search_querry(open('./reads.fasta', 'r'), k_mer, pd.read_csv('./index.dp'))

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

        p = [start_comparison,index_subs]


    return(p)


def seed_and_extend(sequence, querry_found, max_hamming) :
    
    """
        Displays the read ID and its best mapping position and the number of substitions observed
        
        :param sequence: A sequence of nucleotide
        :param querry_found: The output of the function search_querry()
        :param max_hamming: The maximum number of substitutions allowed
        :type sequence: string
        :type querry_found: dictionary
        :type max_hamming: int
        
        
    """
    sequence = sequence[:-1]
    sai_of_kmer = querry_found[0]
    reads = querry_found[-1]


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
                            results = comparison(sequence, read_good_sense, start_comparison, max_hamming, index_kmer)
                            if (results[0]) not in comparison_result.keys(): #If the number of substitutions is not already in the dictionary, create a new key
                                comparison_result[(results[0])] = [(results[1],k_mer_sens)]
                            else : #Else, add the results to the existing key
                                comparison_result[(results[0])].append((results[1],k_mer_sens))
                            list_of_position.append(start_comparison)
            else : # One correspondance only
                if sai_values != '' : 
                    start_comparison = int(sai_values-index_kmer)
                    if start_comparison >= 0 and start_comparison not in list_of_position:
                        results = comparison(sequence, read_good_sense, start_comparison, max_hamming, index_kmer)
                        if (results[0]) not in comparison_result.keys(): #If the number of substitutions is not already in the dictionary, create a new key
                                comparison_result[(results[0])] = [(results[1],k_mer_sens)]
                        else : #Else, add the results to the existing key
                                comparison_result[(results[0])].append((results[1],k_mer_sens))
                        list_of_position.append(start_comparison)

        
        
        #for key in comparison_result.keys():
            #print(len(comparison_result[key][0][0]))
        #print(sorted(comparison_result.items()))
    

        

        """
        change str(value) for key in int
        """
        
#         if bool(comparison_result): #Find the lowest number of substitutions in the dictionary, then find the farthest left alignment
#             #first_value_sort_dic = next(iter(sorted(comparison_result.items(), key=lambda t: len(t[0][0]))))
#             first_value_sort_dic = sorted(comparison_result.items(), key=lambda t: len(t[0][0]))
#             print(first_value_sort_dic)
# 
#             substitution_min = first_value_sort_dic[0]
# 
#             values = first_value_sort_dic[1]
#             min_seq_index = min(values, key = lambda t: t[0])
#             
#             #Print the result
#             if int(substitution_min) <= max_hamming :
#                 print("read " + str(x+1))
#                 print("Nb substitution : " + str(substitution_min) + " - - " + str(min_seq_index))
#         else :
#             print("read " + str(x+1))
#             print("not found")






seed_and_extend(sequence, querry_found, max_hamming)

################TIME COUNT################
end = time.time()

time_to_analyse = end-start
print("TIME FOR ANALYSE : " + str(round(time_to_analyse,3)) + " secondes")


