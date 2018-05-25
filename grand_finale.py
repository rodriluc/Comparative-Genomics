import re
from collections import defaultdict
import math 
#gc, dinucleotide, amino acid


###################################################################
######################### STATISTICS TOOL #########################
###################################################################

def compute_gc():
    sequence = ''
    #ignore = re.compile('^[ \\t]*#.*', re.IGNORECASE)
    with open('GC content Frequency', 'w') as w:
        with open('03.fa.txt', 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line
                    C=sequence.count('C') #count C
                    G=sequence.count('G') #count G
                    total=len(sequence)
                    res=(C+G)/(total) #sum of GC count/total
            w.write('The GC content frequency is '+ str(res))
    
#look for GC content then sum/sequence

def compute_dinucleo():
    sequence = ''
    #dinucleo = ['AG', 'AA', 'AC', 'AT','CG', 'CA', 'CC', 'CT','GG', 'GA', 'GC', 'GT', 'TG', 'TA', 'TC', 'TT']
    #count = 0
    item = 'AG'
    temp_dict = defaultdict(int)
    with open('Dinucleotide Frequency', 'w') as w:
        with open('03.fa.txt', 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line
                    #total = len(sequence)-1
                for item in range(len(sequence)-1):
                    temp_dict[sequence[item:item+2]] +=1
                for k,v in sorted(temp_dict.items()):
                    w.write('Dinucleotide frequency of '+ k+' is '+ str(v)+'\n')
                    
                '''for item in dinucleo: #specify item
                    if item in sequence:
                        count += 1
                        total=len(sequence)-1
                        res=(count)/(total)'''
                        
    #print ('The dinucleotide frequency of AG is ', res, '\n') #do one at a time or print all
# AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT  
       
def compute_aa(): #tanslate ORF list and then comput_aa
    sequence = ''
    list_aa = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W','Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q'] #20aa
    with open('Amino acid Frequency', 'w') as w:
        with open('03.fa.txt.pfa', 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line       
            for item in list_aa:
                i = sequence.count(item)
                res = i/(len(sequence) - 1)
                w.write('Amino acid frequency of '+item+' is ' + str(res) + '\n')
    


##############################################################
######################### ORF FINDER #########################
##############################################################

def complementDNA():
    s = ''
    with open('03.fa.txt', 'r') as f:
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
                    return s #complement
                    
def rev_complementDNA():
    s = ''
    with open('03.fa.txt', 'r') as f:
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
    with open('ORF finder FASTA', 'w') as w:  
        #codon_list = [] #forward strand list
        for i in range(len(complement)-2): #iterates over all possible positions where a codon begin, so all except last 2
            #codon_list.append(complement[i:i+3])
            #codon_list.append(complement)
            #complement.count('ATG')
            #return codon_list
            correct_way = re.compile(r'(?=(ATG(?:...)*?)(?=TAG|TAA|TGA))') #starts with ATG then go to the next stop codon...so on
            count = 0
            for i in correct_way.findall(complement):
                if len(i) > 100: #chose 100 because of Karlin et al. reference
                    #tack Kajetan format help
                    w.write('>ORF_{}\n{}\n'.format(count, ''.join(str(i)))) #took out set so takes gene duplication into account, since findall takes into account overlaps
                    count+=1
                    
            count = 0
            for i in correct_way.findall(reverse_complement):
                if len(i) > 100:       
                    w.write('>ORF_rev_{}\n{}\n'.format(count, ''.join(str(i)))) #write function sooooooooooooo slow don't know why, chnage to print to see results
                    count+=1
            #REGEX FOR THE WIN
             

######################################################################## 
######################### DISTANCE MATRIX TOOL #########################
########################################################################    

def distance_matrix():
#The tool should compute the distance between two genomes from the DNA statistic above. The distance matrix is then used to create a species tree

#matrix one against all


if __name__ == '__main__':
    #print(compute_gc())   
    #print(compute_dinucleo()) 
    #print(compute_aa())
    #print(complementDNA())  
    print(ORF_finder())
    print(distance_matrix())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



