#!/usr/bin/env python2
import re
from collections import defaultdict
import math 
import numpy as np
from string import maketrans
import sys

##############################################################
######################### ORF FINDER #########################
##############################################################

def complementDNA(input_genome):
    s = ''
    with open(input_genome, 'r') as f:
        for line in f:
            for line in f:
                s+=line
        return s #not complement, raw sequence as string
                    
def rev_complementDNA(input_genome):
    s = ''
    with open(input_genome, 'r') as f:
        for line in f:
                if not line.startswith('>'):
                    s += line 
                    for ch in f:
                        if ch=='A':
                            s=s+'T'
                        elif ch=='T':
                            s=s+'A'
                        elif ch=='G':
                            s=s+'C'
                        else:
                            s=s+'G'
                    return s[::-1] #reverse complement                    
    
         
def ORF_finder(input_genome):
#The input should be a genome file in FASTA format; the output should also be a file in FASTA format with separate entries for each ORF gene sequences and unique names identifying these ORFs

    complement = complementDNA(input_genome)
    reverse_complement = rev_complementDNA(input_genome)
    with open(input_genome + '_ORF_finder.fasta', 'w') as w:  
    #stop_codon = ['TAG','TAA','TGA']  
            
        correct_way = re.compile(r'(?=(ATG(?:...)*?)(?:TAG|TAA|TGA))') #starts with ATG then go to the next stop codon...so on
            #for line in correct_way:
            #read backwards solve ^ problem or $
        count = 0
        for i in correct_way.findall(complement):
            if  len(i) > 300: #chose 300 because of Karlin et al. reference
                w.write('>ORF_{}\n{}\n'.format(count, ''.join(str(i)))) #took out set so takes gene duplication into account, since findall takes into account overlaps with correct_way
                count+=1
       
        count = 0
        for i in correct_way.findall(reverse_complement):
            if len(i) > 300: #only ORFs longer then 300 nucleotides 
                w.write('>ORF_rev_{}\n{}\n'.format(count, ''.join(str(i)))) 
                count+=1


###################################################################
######################### STATISTICS TOOL #########################
###################################################################

def compute_gc(input_genome):
    sequence = ''
    with open(input_genome + '_GC_content_Frequency', 'w') as w:
        with open(input_genome, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line
                    C=sequence.count('C') #count C
                    G=sequence.count('G') #count G
                    total=len(sequence)
                    res = float(C+G)/total #sum of GC count/total
            w.write('The GC content frequency is '+ str(res))
            
def compute_nucleo(input_genome):
    sequence = ''
    with open(input_genome + '_Nucleotide_Frequency', 'w') as w:
        with open(input_genome, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line
                    C=sequence.count('C') #count C
                    G=sequence.count('G') #count G
                    A=sequence.count('A') #count A
                    T=sequence.count('T') #count T
                    total=len(sequence)
                    resG = float(G)/total #sum of individual nucleotide/total
                    resC = float(C)/total 
                    resA = float(A)/total
                    resT = float(T)/total
            w.write('The G nucleotide frequency is '+ str(resG)+'\n')
            w.write('The C nucleotide frequency is '+ str(resC)+'\n') 
            w.write('The A nucleotide frequency is '+ str(resA)+'\n') 
            w.write('The T nucleotide frequency is '+ str(resT))  

def compute_dinucleo(input_genome):
    sequence = ''
    temp_dict = defaultdict(int)
    with open(input_genome + '_Dinucleotide_Frequency', 'w') as w:
        with open(input_genome, 'r') as f:
            for line in f:
                undefined=sequence.count('N') 
                if not line.startswith('>'):
                    sequence += line
                for item in range((len(sequence)-1)-(2*undefined)): #excludes N from final count
                    temp_dict[sequence[item:item+2]] +=1 #finds nucleotides*2, puts in dict and counts them
                for k,v in sorted(temp_dict.items()):
                    w.write('Dinucleotide frequency of '+ k+' is '+ str(v)+'\n')

def trans_aa(input_genome): #translate ORF list and then compute_aa
    sequence = ''
    trans_dict = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TGG':'W', 'TGT':'C',
    'TGC':'C'} #didn'include stop since not in ORF

    with open(input_genome, 'r') as f:
        temp_list = []
        for line in f:
            if not line.startswith('>'):
                sequence = line 
                new = ''
                for i in range(0, len(sequence), 3): #looks in window of 3 nucleotides at a time
                    if sequence[i:i+3] in trans_dict.keys():
                        new += trans_dict[sequence[i:i+3]] #translates the found 3 nucleo(key) to the value in trans_dict            
                temp_list.append(new)
        return ''.join(temp_list)
                
                    
def compute_aa(input_genome):
    list_aa = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W','Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q'] 
    trans = trans_aa(input_genome)
    with open(input_fasta + '_amino_acid_frequency', 'w') as w:
        for item in list_aa:
            i = trans.count(item) #counts aa from list_aa found in trans(translated sequence)
            res = float(i)/(len(trans)-1) #number of aa/total
            w.write('Amino acid frequency of '+item+' is ' + str(res) + '\n')


def compute_diaa(input_genome): 
    list_diaa = ['GG', 'GA', 'GL', 'GM', 'GF', 'GW', 'GK', 'GQ', 'GE', 'GS',
        'GP', 'GV', 'GI', 'GC', 'GY', 'GH', 'GR', 'GN', 'GD', 'GT',
        'AG', 'AA', 'AL', 'AM', 'AF', 'AW', 'AK', 'AQ', 'AE', 'AS',
        'AP', 'AV', 'AI', 'AC', 'AY', 'AH', 'AR', 'AN', 'AD', 'AT',
        'LG', 'LA', 'LL', 'LM', 'LF', 'LW', 'LK', 'LQ', 'LE', 'LS',
        'LP', 'LV', 'LI', 'LC', 'LY', 'LH', 'LR', 'LN', 'LD', 'LT',
        'MG', 'MA', 'ML', 'MM', 'MF', 'MW', 'MK', 'MQ', 'ME', 'MS',
        'MP', 'MV', 'MI', 'MC', 'MY', 'MH', 'MR', 'MN', 'MD', 'MT',
        'FG', 'FA', 'FL', 'FM', 'FF', 'FW', 'FK', 'FQ', 'FE', 'FS',
        'FP', 'FV', 'FI', 'FC', 'FY', 'FH', 'FR', 'FN', 'FD', 'FT',
        'WG', 'WA', 'WL', 'WM', 'WF', 'WW', 'WK', 'WQ', 'WE', 'WS',
        'WP', 'WV', 'WI', 'WC', 'WY', 'WH', 'WR', 'WN', 'WD', 'WT',
        'KG', 'KA', 'KL', 'KM', 'KF', 'KW', 'KK', 'KQ', 'KE', 'KS',
        'KP', 'KV', 'KI', 'KC', 'KY', 'KH', 'KR', 'KN', 'KD', 'KT',
        'QG', 'QA', 'QL', 'QM', 'QF', 'QW', 'QK', 'QQ', 'QE', 'QS',
        'QP', 'QV', 'QI', 'QC', 'QY', 'QH', 'QR', 'QN', 'QD', 'QT',
        'EG', 'EA', 'EL', 'EM', 'EF', 'EW', 'EK', 'EQ', 'EE', 'ES',
        'EP', 'EV', 'EI', 'EC', 'EY', 'EH', 'ER', 'EN', 'ED', 'ET',
        'SG', 'SA', 'SL', 'SM', 'SF', 'SW', 'SK', 'SQ', 'SE', 'SS',
        'SP', 'SV', 'SI', 'SC', 'SY', 'SH', 'SR', 'SN', 'SD', 'ST',
        'PG', 'PA', 'PL', 'PM', 'PF', 'PW', 'PK', 'PQ', 'PE', 'PS',
        'PP', 'PV', 'PI', 'PC', 'PY', 'PH', 'PR', 'PN', 'PD', 'PT',
        'VG', 'VA', 'VL', 'VM', 'VF', 'VW', 'VK', 'VQ', 'VE', 'VS',
        'VP', 'VV', 'VI', 'VC', 'VY', 'VH', 'VR', 'VN', 'VD', 'VT',
        'IG', 'IA', 'IL', 'IM', 'IF', 'IW', 'IK', 'IQ', 'IE', 'IS',
        'IP', 'IV', 'II', 'IC', 'IY', 'IH', 'IR', 'IN', 'ID', 'IT',
        'CG', 'CA', 'CL', 'CM', 'CF', 'CW', 'CK', 'CQ', 'CE', 'CS',
        'CP', 'CV', 'CI', 'CC', 'CY', 'CH', 'CR', 'CN', 'CD', 'CT',
        'YG', 'YA', 'YL', 'YM', 'YF', 'YW', 'YK', 'YQ', 'YE', 'YS',
        'YP', 'YV', 'YI', 'YC', 'YY', 'YH', 'YR', 'YN', 'YD', 'YT',
        'HG', 'HA', 'HL', 'HM', 'HF', 'HW', 'HK', 'HQ', 'HE', 'HS',
        'HP', 'HV', 'HI', 'HC', 'HY', 'HH', 'HR', 'HN', 'HD', 'HT',
        'RG', 'RA', 'RL', 'RM', 'RF', 'RW', 'RK', 'RQ', 'RE', 'RS',
        'RP', 'RV', 'RI', 'RC', 'RY', 'RH', 'RR', 'RN', 'RD', 'RT',
        'NG', 'NA', 'NL', 'NM', 'NF', 'NW', 'NK', 'NQ', 'NE', 'NS',
        'NP', 'NV', 'NI', 'NC', 'NY', 'NH', 'NR', 'NN', 'ND', 'NT',
        'DG', 'DA', 'DL', 'DM', 'DF', 'DW', 'DK', 'DQ', 'DE', 'DS',
        'DP', 'DV', 'DI', 'DC', 'DY', 'DH', 'DR', 'DN', 'DD', 'DT',
        'TG', 'TA', 'TL', 'TM', 'TF', 'TW', 'TK', 'TQ', 'TE', 'TS',
'TP', 'TV', 'TI', 'TC', 'TY', 'TH', 'TR', 'TN', 'TD', 'TT']
    trans = trans_aa(input_genome)
    
    with open(input_genome + '_Diamino_acid_Frequency', 'w') as w:
        for item in list_diaa:
            i = trans.count(item) #counts diaa from list_diaa found in trans(translated sequence)
            total_count = len(trans)/6 #divide by 6 because diaa and take into acount that looking at codons before
            res = (float(i)/total_count)*100
            w.write(item+' : '+ str(res)+'\n')            
            
                   
if __name__ == '__main__':
    input_genome = sys.argv[1]
    compute_gc(input_genome)   
    compute_nucleo(input_genome)
    compute_dinucleo(input_genome) 
    compute_diaa(input_genome)
    compute_aa(input_genome) 
    ORF_finder(input_genome)

