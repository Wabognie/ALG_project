"""
ALG Project : mapping without gap, SNPs detection
M2 IBI 2020-2021
Benjamin Blanc
Tiphaine Casy
"""
################IMPORT SECTION################
import pandas as pd
import sys
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
k_mer = 5 ##args.k
max_hamming = 6 ##args.max_hamming

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

        if i_min >=0 and j_max >=0 :
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
#get_querry(get_BWT(sequence)[0], "GG", get_N(sequence),get_BWT(sequence)[-1], sa)

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
    for read_lines in reads :
        read_lines = str(read_lines).replace('\n','')
        if '>' not in read_lines :
            kmer_sai = {}
            k_mer_modified = k_mer

            for x in range(0, len(read_lines)-k_mer_modified+1):
                sai_querry = get_querry(index['BWT'], str(read_lines[x:k_mer_modified]), get_N(sequence), index['SA[i]'])
                if sai_querry != '' and str(read_lines[x:k_mer_modified]) not in kmer_sai.keys() :
                    kmer_sai[str(read_lines[x:k_mer_modified])] = sai_querry ##(sai_querry,x)

                k_mer_modified +=1
            list_querry.append(kmer_sai)


    return(list_querry)
querry_found = search_querry(open('./reads_bis.fasta', 'r'), k_mer, pd.read_csv('./index.dp'))
##search_querry(open(str(args.reads), 'r'), args.k_mer, open(str(args.index), 'r'))


"""
TO DO HERE : read genome with little sa[i] found for querry and compare

BECAREFULL : IF NOT SA[i] FOUND DO REVERSE TRANSCRIPTION OF READ (ex : read 4 and read 10)
"""
def comparison(or_sequence, reads, start_comparison, max_substitution) :
    substitution = 0
    index_subs = []
    comparison_changed = start_comparison
    if start_comparison < len(or_sequence)-len(reads) :
        for x in range(len(reads)-1):
            if or_sequence[comparison_changed] != reads[x]:
                substitution +=1
                result = (or_sequence[start_comparison], start_comparison, reads[x])
                index_subs.append(result)
            comparison_changed+=1

    #print(substitution)
    print(index_subs)
    if substitution == 0 :
        print("Nb of substitution : " + str(substitution) + ' - - ' + str(start_comparison))
        return substitution, start_comparison

    elif substitution < max_substitution :
        print("Nb of substitution : " + str(substitution) + ' - - ' + str(start_comparison))
        return substitution, start_comparison

def seed_and_extend(reads, ref, querry_found, max_hamming) :
    sequence_ref = ''
    for ref_line in ref :
        ref_line = str(ref_line).replace('\n','')
        if '>' not in ref_line :
            sequence_ref = ref_line

    reads_list = []
    for reads_line in reads :
        reads_line = str(reads_line).replace('\n','')
        if '>' not in reads_line :
            reads_list.append(reads_line)

    for x in range(0,len(querry_found)) : ##open each read saved into dictionnary
        #substitution = 0
        itera = next(iter(querry_found[x])) ##take the first value of list so the first k_mer, it's can maybe an other
        #print(querry_found[x])
        #print(itera)
        #print(querry_found[x][itera])
        #### test sur le premier kmer trouve
        for i in range(len(querry_found[x][itera])) :
            comparison(str(sequence_ref), str(reads_list[x]),querry_found[x][itera][i], max_hamming)






        #for key in querry_found[x].keys() : ##open each sa[i] values of each k_mer for each reads
            #print("KEY : " + str(key) + " ** " + str(querry_found[x][key]))
            #print("BLABLABLA : " + str(next(iter(querry_found[x]))))
        #print(len(querry_found[x]))

seed_and_extend(open('./reads_bis.fasta', 'r'), open('./reference.fasta','r'),querry_found, max_hamming)

################TIME COUNT################
end = time.time()

time_to_analyse = end-start
print("TIME FOR ANALYSE : " + str(round(time_to_analyse,3)) + " secondes")



"""

ADD NUMBER OF K_MER TO KNOW WHERE WE ARE INTO THE ORIGIN sequence
SUBSTRACT NB - start_comparison
"""





"""
Recuperation des reads dans le fichiers, etrecuperation de la table de burrow wheller (faite avec le index.py)
k-mer : taille definie par l'utilisateur
prendre tous les k-mer du read, le chercher dans la BWT en la remontant et return le sa[i] correspondant (PT : OPTIMISER)
    dictionnaire : nb_kmer : [SA[i]]
    attention sa[i]<taille du genome-taille du read (condition a respecter)
 pour tous les indexs du k-mer etend la lecture avec une comparaison de tous les nucleotides un a un, forward and reverse
    si difference entre read et genome ref = substitution
comptage des substitutions
    choix de l'indexe avec le moins de substitutions, s'il en existe plusieurs c'est la position la plus a gauche dans le genome et donc l'indexe le plus faible

    stockage des substitutions qqpart pour les reutiliser pour le document final

retourner l'alignement entre les deux : entre genome et reads
stockage des infos dans le fichier vcf
"""
