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


################OPEN READS FILE################
file = open('./reads.fasta', 'r')
k_mer = 5



################BWT INFORMATION AND RESEARCH################
def get_N(sequence):
    N = {}
    for i in sequence :
        if i not in N.keys():
            N[i]=1
        else :
            N[i]+=1
    return N
#get_N("GGCGGCACCGC$")

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
#LF("G",3,get_N(sequence))

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
#R_table(4,get_BWT(sequence)[0])

def get_S(BWT,N):
    i = 0
    S = BWT[i]+"$"

    for i in range(len(BWT)-2):
        next_i = LF(S[0],R_table(int(i),BWT), N)
        S = BWT[next_i] + S
        i = next_i
        return S
#get_S(get_BWT(sequence)[0], get_N(sequence))

def get_querry(BWT, Q, N, BWT_list, sa) :
    last_char = Q[-1]
    i = LF(last_char,1,N)
    j = LF(last_char, N[last_char], N)
    i_min = 0
    j_max = 0
    for x in reversed(range(len(Q)-1)):
        index = []
        for l in range(i, j+1):
            if str(BWT[l]) == str(Q[x]) :
                index.append(l)
                i_min = min(index)
                j_max = max(index)

        i = LF(BWT[i_min],R_table(i_min, BWT),N)
        j = LF(BWT[j_max],R_table(j_max, BWT),N)

        print(i)
        print(j)
        if int(i) == int(j) :
            print(sa[i])
        else :
            print("SA I : " + str(sa[i]) + "/ / "+ str(sa[j]))
#get_querry(get_BWT(sequence)[0], "GG", get_N(sequence),get_BWT(sequence)[-1], sa)

################OPEN READS FILE################
for line in file :
    line = str(line).replace('\n','')
    if '>' not in line :
        print(line)
        for x in range(0, len(line)-k_mer+1):
            print(line[x:k_mer])
            #get_querry(get_BWT(sequence)[0], str(line[x:k_mer]), get_N(sequence),get_BWT(sequence)[-1], sa)
            k_mer +=1


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
