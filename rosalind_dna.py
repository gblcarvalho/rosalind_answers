import re
import urllib.request
import string
import itertools

# DNA
def counting_dna_nucleotides(dna):
    return (dna.count("A"), dna.count("C"), dna.count("G"), dna.count("T"))

# RNA
def transcribing_dna_into_rna(dna):
    return dna.replace("T","U")

# REVC
def complementing_a_string_of_dna(dna):
    return "".join([ {"A":"T", "T":"A", "C":"G", "G":"C"}[n] for n in dna])[::-1]

# IPRB
def mendels_first_law(k,m,n):
    total = k + m + n
    Kp = k/total
    Mp = ((m/total) * (k/(total-1))) + ((m/total) * ((m-1)/(total-1)) * 0.75) + ((m/total) * (n/(total-1)) * 0.5)
    Np = ((n/total) * (k/(total-1))) + ((n/total) * (m/(total-1)) * 0.5)
    return Kp + Mp + Np

# FIB
def rabbits_and_recurrense_relations(n , k):
    i = 1
    C = 0
    c = 0
    while i <= n:
        if i == 1:
            C = 0
            c = 1
        else:
            C_aux = C
            C = C + c
            c = C_aux * k
        i += 1
    return C + c

# GC
def read_fasta(fasta):
    dnas = {}
    arquivo = fasta.split(">")
    for linha in arquivo:
        ls = linha.strip()
        if ls:
            l = ls.split('\n')
            label = l[0]
            dna = "".join(l[1:])
            dnas[label] = dna
    return dnas

def gc_content(dna):
    return ((dna.count("G") + dna.count("C")) * 100) / len(dna)

def max_gc_content(filename):
    dnas = read_fasta(filename)
    max_label = ""
    max_dna = ""
    max_gc = 0
    for label, dna in dnas.items():
        gc = gc_content(dna)
        if gc > max_gc:
            max_gc = gc
            max_label = label
            max_dna = dna
    return max_label, max_dna, max_gc

# PROT
def translating_rna_into_protein(rna):
    table = {'UUU':'F', 'CUU':'L', 'AUU':'I', 
             'GUU':'V', 'UUC':'F', 'CUC':'L', 
             'AUC':'I', 'GUC':'V', 'UUA':'L', 
             'CUA':'L', 'AUA':'I', 'GUA':'V', 
             'UUG':'L', 'CUG':'L', 'AUG':'M', 
             'GUG':'V', 'UCU':'S', 'CCU':'P', 
             'ACU':'T', 'GCU':'A', 'UCC':'S', 
             'CCC':'P', 'ACC':'T', 'GCC':'A', 
             'UCA':'S', 'CCA':'P', 'ACA':'T', 
             'GCA':'A', 'UCG':'S', 'CCG':'P', 
             'ACG':'T', 'GCG':'A', 'UAU':'Y', 
             'CAU':'H', 'AAU':'N', 'GAU':'D', 
             'UAC':'Y', 'CAC':'H', 'AAC':'N', 
             'GAC':'D', 'UAA':'Stop', 'CAA':'Q', 
             'AAA':'K', 'GAA':'E', 'UAG':'Stop', 
             'CAG':'Q', 'AAG':'K', 'GAG':'E', 
             'UGU':'C', 'CGU':'R', 'AGU':'S', 
             'GGU':'G', 'UGC':'C', 'CGC':'R', 
             'AGC':'S', 'GGC':'G', 'UGA':'Stop', 
             'CGA':'R', 'AGA':'R', 'GGA':'G', 
             'UGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G'}
    regex = re.compile(r'.{3}')
    protein = ""
    for n in regex.finditer(rna):
        p = table[n.group()]
        if p == "Stop":
            break
        protein += p
    return protein

# SUBS
def finding_a_motif_in_dna(motif, dna):
    regex = re.compile(r'(?=%s)' % motif)
    return [ m.start() + 1 for m in regex.finditer(dna)]

#HAMM
def counting_point_mutations(dna1, dna2):
    i = 0
    total = 0
    while i < len(dna1):
        if dna1[i] != dna2[i]:
            total += 1
        i += 1
    return total

#GRAPH
def overlap_graphs(dnas, k):
    adj_list = []
    for label, dna in dnas.items():
        for label2, dna2 in dnas.items():
            if label != label2:
                if dna[-k:] == dna2[:k]:
                    adj_list.append([label, label2])
    return adj_list

#MPRT
def fetch_uniprot_fasta(protein):
    fasta = urllib.request.urlopen("http://www.uniprot.org/uniprot/" + protein + ".fasta")
    return fasta.read().decode("utf-8")

def finding_a_protein_motif(protein):
    regex = re.compile(r'(?=N[^P][ST][^P])')
    return [m.start() + 1 for m in regex.finditer(protein)]

def show_a_protein_motif_position(proteins):
    for p in proteins:
        dnas = read_fasta(fetch_uniprot_fasta(p))
        for label, dna in dnas.items():
            positions = finding_a_protein_motif(dna)
            if positions:
                print(p)
                print(" ".join([ str(n) for n in positions ]))

#CONS
def consensus_and_profile(dnas):
    profile = {} 
    for p in ["A", "C", "T", "G"]:
        profile[p] = [0 for i in dnas[0]]

    for dna in dnas:
        for i, n in enumerate(dna):
            profile[n][i] += 1

    i = 0
    consensus = ""
    while i < len(profile["A"]):
        max_value = 0
        p = ""
        if profile["A"][i] > max_value:
            p = "A"
            max_value = profile["A"][i]
        if profile["T"][i] > max_value:
            p = "T"
            max_value = profile["T"][i]
        if profile["C"][i] > max_value:
            p = "C"
            max_value = profile["C"][i]
        if profile["G"][i] > max_value:
            p = "G"
            max_value = profile["G"][i]
        consensus += p
        i += 1
    return profile, consensus

#SPLC
def delete_instrons(dna, introns):
    points = []
    for i in introns:
        regex = re.compile(r'%s' % i)
        for m in regex.finditer(dna):
            for n in range(m.start(), m.start() + len(i)):
                points.append(n)
    dna_res = ""
    for i, l in enumerate(dna):
        if i not in points:
            dna_res += l
    return dna_res

def read_fasta_list(fasta):
    dnas = []
    arquivo = fasta.split(">")
    for linha in arquivo:
        ls = linha.strip()
        if ls:
            l = ls.split('\n')
            dna = "".join(l[1:])
            dnas.append(dna)
    return dnas

#FIBD
def mortal_fibonacci_rabbits(n , m):
    i = 1
    C = 0
    c = [0]
    while i <= n:
        if i == 1:
            C = 0
            c.append(1)
        else:
            C_aux = C
            if i <= m:
                C = C + c[i-1]
            else:
                C = C + c[i-1] - (c[i-m])
            c.append(C_aux)
        i += 1
    return C + c[i-1]

#LCSM
def longest_common_substring(string1, string2):
    longest_substring = ""
    initial = 0
    length = 1
    max_length = len(string1) 
    while True:
        if string1[initial:initial+length] in string2:
            longest_substring = string1[initial:initial+length]
            length += 1
        else:
            initial += 1
            length = len(longest_substring)
        if initial + length > max_length:
            break
    return longest_substring


def finding_a_shared_motif(dnas):
    current_lcs = ""
    total_lcs = set()
    for i in range(1, len(dnas)):
        current_lcs = longest_common_substring(dnas[0], dnas[i])
        for j in range(1, len(dnas)):
            if not current_lcs in dnas[j]:
                current_lcs = longest_common_substring(current_lcs, dnas[j])
        total_lcs.add(current_lcs)

    return max(total_lcs, key=len)

#SSEQ
def subsequence_indexes(palavra, sub):
    indexes = []
    init = 0
    for l in sub:
        init = palavra.find(l,init) + 1

        indexes.append(init)
    return indexes

#LCSQ
def montar_matriz_lcs(palavra1, palavra2):

    matriz = [[ 0 for i in range(len(palavra2) + 1)] for j in range(len(palavra1) + 1)]
    
    for i, a in enumerate(palavra1):
        for j, b in enumerate(palavra2):
            if a == b:
                matriz[i+1][j+1] = matriz[i][j] + 1
            else:
                matriz[i+1][j+1] = max(matriz[i][j+1], matriz[i+1][j])

    return matriz


def descobrir_sequencia(matriz, palavra1, palavra2):
    i = len(palavra1) - 1
    j = len(palavra2) - 1
    sec = ""
    while i >= 0 and j >= 0:
        if palavra1[i] == palavra2[j]:
            sec = sec + palavra1[i]
            i -= 1
            j -= 1
        elif matriz[i][j+1] >= matriz[i+1][j]:
            i -= 1
        else:
            j -= 1

    return sec[::-1]
    

def lcs(palavra1, palavra2):
    matriz = montar_matriz_lcs(palavra1, palavra2)
    return descobrir_sequencia(matriz, palavra1, palavra2)

def finding_a_shared_spliced_motif(palavra1, palavra2):
    return lcs(palavra1, palavra2)

#ORF
def open_reading_frames(rna):
    table = {'UUU':'F', 'CUU':'L', 'AUU':'I', 
             'GUU':'V', 'UUC':'F', 'CUC':'L', 
             'AUC':'I', 'GUC':'V', 'UUA':'L', 
             'CUA':'L', 'AUA':'I', 'GUA':'V', 
             'UUG':'L', 'CUG':'L', 'AUG':'M', 
             'GUG':'V', 'UCU':'S', 'CCU':'P', 
             'ACU':'T', 'GCU':'A', 'UCC':'S', 
             'CCC':'P', 'ACC':'T', 'GCC':'A', 
             'UCA':'S', 'CCA':'P', 'ACA':'T', 
             'GCA':'A', 'UCG':'S', 'CCG':'P', 
             'ACG':'T', 'GCG':'A', 'UAU':'Y', 
             'CAU':'H', 'AAU':'N', 'GAU':'D', 
             'UAC':'Y', 'CAC':'H', 'AAC':'N', 
             'GAC':'D', 'UAA':'Stop', 'CAA':'Q', 
             'AAA':'K', 'GAA':'E', 'UAG':'Stop', 
             'CAG':'Q', 'AAG':'K', 'GAG':'E', 
             'UGU':'C', 'CGU':'R', 'AGU':'S', 
             'GGU':'G', 'UGC':'C', 'CGC':'R', 
             'AGC':'S', 'GGC':'G', 'UGA':'Stop', 
             'CGA':'R', 'AGA':'R', 'GGA':'G', 
             'UGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G'}
    
    start_condons = [i.start() for i in re.finditer('AUG', rna)]

    regex = re.compile(r'.{3}')
    proteins = set()


    for sc in start_condons:
        protein = ""
        for n in regex.finditer(rna[sc:]):
            p = table[n.group()]
            if p == "Stop":
                if protein:
                    proteins.add(protein)
                break
            protein += p
    return proteins

def distinct_protein_candidates(dna):
    rna = transcribing_dna_into_rna(dna)
    dna = complementing_a_string_of_dna(dna)
    rna2 = transcribing_dna_into_rna(dna)

    return open_reading_frames(rna) | open_reading_frames(rna2)

#PERM
def permutations(number):
    return [p for p in itertools.permutations(range(1,number+1))]


#REVP
def reverse_palindrome(dna):
    min_length = 4
    max_length = 12
    l = min_length
    i = 0
    positions = []
    while True:
        if i + min_length > len(dna):
            break
        if dna[i:i+l] == complementing_a_string_of_dna(dna[i:i+l]):
            positions.append((i+1,l))
        if l == max_length or i + l >= len(dna):
            l = min_length
            i += 1
        else:
            l += 1
    return positions

#PRTM
def calculating_protein_mass(protein):
    table =  {'A': 71.03711,
                'C': 103.00919,
                'D': 115.02694,
                'E': 129.04259,
                'F': 147.06841,
                'G': 57.02146,
                'H': 137.05891,
                'I': 113.08406,
                'K': 128.09496,
                'L': 113.08406,
                'M': 131.04049,
                'N': 114.04293,
                'P': 97.05276,
                'Q': 128.05858,
                'R': 156.10111,
                'S': 87.03203,
                'T': 101.04768,
                'V': 99.06841,
                'W': 186.07931,
                'Y': 163.06333 }

    return sum([table[a] for a in protein])

    


        
