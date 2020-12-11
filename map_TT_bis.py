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
k_mer = 19 ##args.k
max_hamming = 5 ##args.max_hamming

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
    Function used to research all sa[i] k_mer position into reference genome
    return a list of list with different elements
    """
    list_querry = []
    read_information = {}
    for read_lines in reads :
        read_lines = str(read_lines).replace('\n','')
        reverse_read_lines = reverse_transcript(str(read_lines))

        if '>' not in read_lines :
            kmer_sai = {}
            k_mer_modified = k_mer
            if "+" not in read_information.keys() and "-" not in read_information.keys():
                read_information['+'] = [str(read_lines)]
                read_information['-'] = [str(reverse_read_lines)]
            else :
                read_information['+'].append(str(read_lines))
                read_information['-'].append(str(reverse_read_lines))



            for x in range(0, len(read_lines)-k_mer_modified+1):
                sai_querry = get_querry(index['BWT'], str(read_lines[x:k_mer_modified]), get_N(sequence), index['SA[i]'])

                if sai_querry != '':
                    kmer_sai[str(read_lines[x:k_mer_modified]),x] = [sai_querry,x, "+"] ##(sai_querry,x)

                else :
                    sai_querry = get_querry(index['BWT'],str(reverse_read_lines[x:k_mer_modified]), get_N(sequence), index['SA[i]'])
                    kmer_sai[str(reverse_read_lines[x:k_mer_modified]),x] = [sai_querry,x, "-"] ##(sai_querry,x)
                    #if sai_querry == '':
                        #kmer_sai[str(reverse_read_lines[x:k_mer_modified]),x] = ["NONE",x, "None"]

                k_mer_modified +=1
            list_querry.append(kmer_sai)

    return(list_querry,read_information)
querry_found = search_querry(open('./reads_bis.fasta', 'r'), k_mer, pd.read_csv('./index.dp'))

def comparison(or_sequence, reads, start_comparison, max_substitution, index_k_mer) :
    substitution = max_substitution
    index_subs = []

    comparison_changed = start_comparison
    if start_comparison <= len(or_sequence)-len(reads) :
        substitution = 0
        for x in range(len(reads)-1):
            if or_sequence[comparison_changed] != reads[x]:
                substitution +=1
                result = (or_sequence[start_comparison], start_comparison, reads[x])
                index_subs.append(result)
            comparison_changed+=1
        p = [substitution, start_comparison]


    return(p)


def seed_and_extend(sequence, querry_found, max_hamming) :
    sequence = sequence[:-1]
    sai_of_kmer = querry_found[0]
    reads = querry_found[-1]


    for x in range(0,len(sai_of_kmer)) : ##pour tous les reads
        comparison_result = {}
        for y in range(0,len(sai_of_kmer[x])): ##pour tous les k_mer de chaque reads
            sai_values = list(sai_of_kmer[x].values())[y][0]
            index_kmer = list(sai_of_kmer[x].values())[y][1]
            k_mer_sens = list(sai_of_kmer[x].values())[y][-1]

            read_good_sense = reads[k_mer_sens][x]
            if isinstance(sai_values, list): ##verify if sai is a list or not
                for i in sai_values :
                    if i != '':
                        start_comparison = int(i-index_kmer)
                        if start_comparison >= 0 :
                            results = comparison(sequence, read_good_sense, start_comparison, max_hamming, index_kmer)
                            if str(results[0]) not in comparison_result.keys():
                                comparison_result[str(results[0])] = [(int(results[1]),k_mer_sens)]
                            else :
                                comparison_result[str(results[0])].append((int(results[1]),k_mer_sens))
            else :
                if sai_values != '' :
                    start_comparison = int(sai_values-index_kmer)
                    if start_comparison >= 0 :
                        results = comparison(sequence, read_good_sense, start_comparison, max_hamming, index_kmer)
                        if str(results[0]) not in comparison_result.keys():
                            comparison_result[str(results[0])] = [(int(results[1]),k_mer_sens)]
                        else :
                            comparison_result[str(results[0])].append((int(results[1]),k_mer_sens))


        #print(comparison_result)
        #print(sorted(comparison_result.items()))

        """
        change str(value) for key in int
        """
        first_value_sort_dic = next(iter(sorted(comparison_result.items())))

        substitution_min = first_value_sort_dic[0]

        values = first_value_sort_dic[1]
        min_seq_index = min(values, key = lambda t: t[0])

        print("Nb substitution : " + str(substitution_min) + " - - " + str(min_seq_index))






seed_and_extend(sequence, querry_found, max_hamming)

################TIME COUNT################
end = time.time()

time_to_analyse = end-start
print("TIME FOR ANALYSE : " + str(round(time_to_analyse,3)) + " secondes")
