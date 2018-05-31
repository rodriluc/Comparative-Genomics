import re
from collections import defaultdict
import math 
import numpy as np
from Bio.Data import CodonTable #not in use
import itertools
from collections import Counter
#gc, dinucleotide, amino acid


###################################################################
######################### STATISTICS TOOL #########################
###################################################################

def compute_gc():
    sequence = ''
    #ignore = re.compile('^[ \\t]*#.*', re.IGNORECASE)
    with open('GC content Frequency: 50.fa.txt', 'w') as w:
        with open('50.fa.txt', 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line                
                    
                    C=sequence.count('C') #count C
                    G=sequence.count('G') #count G
                    total=len(sequence)
                    #print(total)
                    res = float(C+G)/total #sum of GC count/total
                 #print(res)
            w.write('The GC content frequency is '+ str(res))
            #return res
#look for GC content then sum/sequence

def compute_nucleo():
    sequence = ''
    with open('Nucleotide Frequency: 50.fa.txt', 'w') as w:
        with open('50.fa.txt', 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line
                    C=sequence.count('C') 
                    G=sequence.count('G') 
                    A=sequence.count('A') 
                    T=sequence.count('T')
                    total=len(sequence)
                    resG = float(G)/total
                    resC = float(C)/total 
                    resA = float(A)/total
                    resT = float(T)/total
            w.write('The G nucleotide frequency is '+ str(resG)+'\n')
            w.write('The C nucleotide frequency is '+ str(resC)+'\n') 
            w.write('The A nucleotide frequency is '+ str(resA)+'\n') 
            w.write('The T nucleotide frequency is '+ str(resT))     

def compute_dinucleo():
    sequence = ''
    #dinucleo = ['AG', 'AA', 'AC', 'AT','CG', 'CA', 'CC', 'CT','GG', 'GA', 'GC', 'GT', 'TG', 'TA', 'TC', 'TT']
    #count = 0
    #item = 'AG'
    temp_dict = defaultdict(int)
    with open('Dinucleotide Frequency: 03.fa.txt', 'w') as w:
        with open('03.fa.txt', 'r') as f:
            for line in f:
                undefined=sequence.count('N')
                if not line.startswith('>'):
                    sequence += line
                    
                    #total = len(sequence)-1
                for item in range((len(sequence)-1)-(2*undefined)):
                    temp_dict[sequence[item:item+2]] +=1
                for k,v in sorted(temp_dict.items()):
                    total = sum(temp_dict.values())
                    result = float(v)/total
                    #print(result)
                    w.write('Dinucleotide frequency of '+ k+' is '+ str(result)+'\n')
                    
                '''for item in dinucleo: #specify item
                    if item in sequence:
                        count += 1
                        total=len(sequence)-1
                        res=(count)/(total)'''
                        
    #print ('The dinucleotide frequency of AG is ', res, '\n') #do one at a time or print all
# AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT  
       
def trans_aa(): #tanslate ORF list and then comput_aa
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
    #table = CodonTable.ambiguous_dna_by_id[1]

    with open('ORF finder FASTA: 50.fa.txt', 'r') as f:
        temp_list = []
        for line in f:
            if not line.startswith('>'):
                sequence = line 
                new = ''
                for i in range(0, len(sequence), 3):
                    #print(sequence[i:i+3])
                    if sequence[i:i+3] in trans_dict.keys():
                        new += trans_dict[sequence[i:i+3]] #if have to create own dictionary              
                        #new_list.extend(new)    
                temp_list.append(new)
        return ''.join(temp_list)
                
                    #protein_seq = _translate_str(sequence, table)
                    
def convert_protein():

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
    'TGC':'C'} 
    
    with open('ORF finder FASTA: 03.fa.txt', 'r') as f:
        with open('proteome_03', 'w') as w:
        
            '''prot_dict = {}
            names = []
            temp_list = []
            for line in f:
                #print(line)
                if line.startswith('>'):
                    names.append(line[0:-1])
                    #print(names)
                if not line.startswith('>'):
                    sequence = line 
                    new = ''
                    for i in range(0, len(sequence), 3):
                        #print(sequence[i:i+3])
                        if sequence[i:i+3] in trans_dict.keys():
                            new += trans_dict[sequence[i:i+3]]     
                    #temp_list.append(''.join(new))
                    w.write(names+'\n'+new) '''
                    #return prot_dict
            for line[1::2] in f:
                sequence = line 
                new = ''
                for i in range(0, len(sequence), 3):
                    #print(sequence[i:i+3])
                    if sequence[i:i+3] in trans_dict.keys():
                        new += trans_dict[sequence[i:i+3]] 
                w.write(line[0::2]+'\n'+new)
                

                    
def compute_aa():
    list_aa = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W','Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q'] 
    trans = trans_aa()
    #print(trans)
    with open('Amino acid Frequency: 03.fa.txt', 'w') as w:
        for item in list_aa:
            i = trans.count(item)
            #print (len(trans), i, i/len(trans))
            res = float(i)/(len(trans)-1)
            #print(res)
            w.write('Amino acid frequency of '+item+' is ' + str(res) + '\n')
    
def compute_diaa(): 
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
    trans = trans_aa()
    #print(trans)
    
    with open('Diamino acid Frequency: 50.fa.txt', 'w') as w:
        for item in list_diaa:
            i = trans.count(item)
            #print (len(trans), i, i/len(trans))
            res = float(i)/(len(trans)-1) #does not count the X
            #print(res)
            w.write(item+' : '+ str(res)+'\n')
            
    '''temp_dict = defaultdict(int)
    with open('Diamino acid Frequency: 03.fa.txt', 'w') as w:
        for line in trans:   
            for item in range(len(trans)):
                temp_dict[trans[item:item+2]] +=1
                for k,v in sorted(temp_dict.items()):
                    total = sum(temp_dict.values())
                    result = float(v)/total
                    #print(result)
                    w.write(k+' : '+ str(result)+'\n')'''
                    
                    
    #option 1 works BUT very slowly, test option 2
                        
    '''list_aa1 = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W','Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q'] 
    list_aa2 = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W','Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q']
    zip_diaa = [zip(x,list_aa2) for x in itertools.permutations(list_aa1,len(list_aa2))]
    #print (zip_diaa)
    #diaa = [str(x)+str(y) for x in list_aa1 for y in list_aa2]
    #print (diaa)
    poss_combo = dict(enumerate(zip_diaa,1)) 
    kl = []
    templ = []
    with open('Diamino acid Frequency: 03.fa.txt', 'w') as w:
        for k, v in zip_diaa.items():
            kl.append(v)
        set_keys = set(kl)
        for item in range(0, len(trans)):
            read_din = trans[i:i+2]
            templ.append(read_din)
        #print(kl)
        diaa_dict = Counter(templ)
        key_list = (set_keys).intersection(diaa_dict)
        for k,v in diaa_dict.items():
            for item in key_list:
                if k==item:
                    result = float(v/(len(diaa_dict))
                #w.write(k + ' : ' + str(result)+'\n')
                    k, result'''        
       

##############################################################
######################### ORF FINDER #########################
##############################################################

def complementDNA():
    s = ''
    with open('50.fa.txt', 'r') as f:
        for line in f:
            s+=line
        return s #complement USELESS'''
                    
def rev_complementDNA():
    s = ''
    with open('50.fa.txt', 'r') as f:
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
    
def ORF_finder():
#The input should be a genome file in FASTA format; the output should also be a file in FASTA format with separate entries for each ORF gene sequences and unique names identifying these ORFs
    
    complement = complementDNA()
    reverse_complement = rev_complementDNA()
    #stop_codon = ['TAG','TAA','TGA']
    with open('ORF finder FASTA: 50.fa.txt', 'w') as w:  
        #codon_list = [] #forward strand list
        #for i in range(len(complement)): #iterates over all possible positions where a codon begin, so all except last 2
            #print(len(complement))
            #codon_list.append(complement[i:i+3]) 
            #codon_list.append(complement)
            #complement.count('ATG')
            #return codon_list
            
        correct_way = re.compile(r'(?=(ATG(?:...)*?)(?:TAG|TAA|TGA))') #starts with ATG then go to the next stop codon...so on
            #for line in correct_way:
                
            #return set(correct_way.findall(complement))
            #i need to exclude the strings with start codons inside
            #read backwards solve ^ problem or $
        count = 0
        for i in correct_way.findall(complement):
            if  len(i) > 300: #chose 100 because of Karlin et al. reference
                #tack Kajetan format help
                w.write('>ORF_{}\n{}\n'.format(count, ''.join(str(i)))) #took out set so takes gene duplication into account, since findall takes into account overlaps...way too many kept the set
                count+=1
        
        #break        
        count = 0
        for i in correct_way.findall(reverse_complement):
            if len(i) > 300:        #100<= len(i) <=500:
                w.write('>ORF_rev_{}\n{}\n'.format(count, ''.join(str(i)))) #write function sooooooooooooo slow don't know why, change to print to see results
                count+=1
                
            #break
            #REGEX FOR THE WIN
            #print(string, file=filename)

######################################################################## 
######################### DISTANCE MATRIX TOOL #########################
########################################################################    

def distance_matrix_gc():
#The tool should compute the distance between two genomes from the DNA statistic above. The distance matrix is then used to create a species tree

#matrix one against all
# numpy array of 2 vectors
    #GC = compute_gc()
    '''path = ~/Documents/Comparative\Genomics/Final\Project
    for genomes in path
        if genomes.endswith('.fa.txt'):
            genome_files.append(genomes)
            with open (path + genome_files) as f:
                for line in f:
                    if not line.startswith('>'):
                        sequence += line
                        C=sequence.count('C') 
                        G=sequence.count('G') 
                        total=len(sequence)
                        res = float(C+G)/total'''
    sequence=''  
    with open('Distance Matrix', 'w') as w:                  
        with open ('03.fa.txt') as f1:
            with open ('28.fa.txt') as f2:
                with open ('43.fa.txt') as f3:
                    with open ('48.fa.txt') as f4:
                        with open ('50.fa.txt') as f5: 
                            for line in f1:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res1 = float(CG)/total
                            for line in f2:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res2 = float(CG)/total
                            for line in f3:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res3 = float(CG)/total
                            for line in f4:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res4 = float(CG)/total
                            for line in f5:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res5 = float(CG)/total
        #print (res1)
        d_matrix = np.zeros((5,5))
        #print(d_matrix)
        input_genomes = [res1, res2, res3, res4, res5]
        for x in range(0,len(input_genomes)):
            for y in range(0,len(input_genomes)):
                d_matrix[x,y] = math.sqrt((input_genomes[x]-input_genomes[y])**2)
    print(d_matrix)
    
def distance_matrix_dinucl():
    sequence=''  
    with open('Distance Matrix', 'w') as w:                  
        with open ('03.fa.txt') as f1:
            with open ('28.fa.txt') as f2:
                with open ('43.fa.txt') as f3:
                    with open ('48.fa.txt') as f4:
                        with open ('50.fa.txt') as f5:  
                            temp_dict = defaultdict(int)
                            for line in f1:
                                if not line.startswith('>'):
                                    sequence += line
                                    for item in range(len(sequence)-1):
                                        temp_dict[sequence[item:item+2]] +=1
                                    for k,v in sorted(temp_dict.items()):
                                        total = sum(temp_dict.values())
                                        res1 = float(v)/total
                            for line in f2:
                                if not line.startswith('>'):
                                    sequence += line
                                    for item in range(len(sequence)-1):
                                        temp_dict[sequence[item:item+2]] +=1
                                    for k,v in sorted(temp_dict.items()):
                                        total = sum(temp_dict.values())
                                        res2 = float(v)/total
                            for line in f3:
                                if not line.startswith('>'):
                                    sequence += line
                                    for item in range(len(sequence)-1):
                                        temp_dict[sequence[item:item+2]] +=1
                                    for k,v in sorted(temp_dict.items()):
                                        total = sum(temp_dict.values())
                                        res3 = float(v)/total
                            for line in f4:
                                if not line.startswith('>'):
                                    sequence += line
                                    for item in range(len(sequence)-1):
                                        temp_dict[sequence[item:item+2]] +=1
                                    for k,v in sorted(temp_dict.items()):
                                        total = sum(temp_dict.values())
                                        res4 = float(v)/total
                            for line in f5:
                                if not line.startswith('>'):
                                    sequence += line
                                    for item in range(len(sequence)-1):
                                        temp_dict[sequence[item:item+2]] +=1
                                    for k,v in sorted(temp_dict.items()):
                                        total = sum(temp_dict.values())
                                        res5 = float(v)/total
        #print (res1)
        d_matrix = np.zeros((5,5))
        #print(d_matrix)
        input_genomes = [res1, res2, res3, res4, res5]
        for x in range(0,len(input_genomes)):
            for y in range(0,len(input_genomes)):
                d_matrix[x,y] = math.sqrt((input_genomes[x]-input_genomes[y])**2)
    print(d_matrix)
    
    
    
if __name__ == '__main__':
    #print(compute_gc())   
    #print(compute_dinucleo()) 
    #print(compute_nucleo())
    #print(trans_aa())
    #print(convert_protein())
    #print(compute_aa())
    print(compute_diaa())
    #print(complementDNA())  
    #print(ORF_finder())
    #print(distance_matrix_gc())
    #print(distance_matrix_dinucl())
    

    
    
    
    
    



