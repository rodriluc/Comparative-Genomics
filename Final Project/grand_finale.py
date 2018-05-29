import re
from collections import defaultdict
import math 
import numpy as np
from string import maketrans
#gc, dinucleotide, amino acid


###################################################################
######################### STATISTICS TOOL #########################
###################################################################

def compute_gc():
    sequence = ''
    #ignore = re.compile('^[ \\t]*#.*', re.IGNORECASE)
    with open('GC content Frequency: 50.fa.txt', 'w') as w:
        with open('50.fa.txt', 'r') as f:
    #with open(filename, 'r') as f:
            #f = f.readlines()
            for line in f:
                #print(line)
                if not line.startswith('>'):
                    sequence += line
                    #print(sequence)
                    C=sequence.count('C') #count C
                    G=sequence.count('G') #count G
                    total=len(sequence)
                    #print(total)
                    res = float(C+G)/total #sum of GC count/total
                 #print(res)
            w.write('The GC content frequency is '+ str(res))
            #return res
#look for GC content then sum/sequence

def compute_dinucleo():
    sequence = ''
    #dinucleo = ['AG', 'AA', 'AC', 'AT','CG', 'CA', 'CC', 'CT','GG', 'GA', 'GC', 'GT', 'TG', 'TA', 'TC', 'TT']
    #count = 0
    item = 'AG'
    temp_dict = defaultdict(int)
    with open('Dinucleotide Frequency: 50.fa.txt', 'w') as w:
        with open('50.fa.txt', 'r') as f:
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
        with open('ORF Finder FASTA', 'r') as f:
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
                    return s #complement
                    
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
        for i in range(len(complement)-2): #iterates over all possible positions where a codon begin, so all except last 2
            #codon_list.append(complement[i:i+3])
            #codon_list.append(complement)
            #complement.count('ATG')
            #return codon_list
            
            correct_way = re.compile(r'(?=(ATG(?:...)*?)(?:TAG|TAA|TGA))') #starts with ATG then go to the next stop codon...so on
            
            #return set(correct_way.findall(complement))
            #i need to exclude the strings with start codons inside
            #read backwards solve ^ problem or $
            count = 0
            for i in correct_way.findall(complement):
                if len(i) > 100: #chose 100 because of Karlin et al. reference
                    #tack Kajetan format help
                    w.write('>ORF_{}\n{}\n'.format(count, ''.join(str(i)))) #took out set so takes gene duplication into account, since findall takes into account overlaps...way too many kept the set
                    count+=1
            
            #break        
            count = 0
            for i in correct_way.findall(reverse_complement):
                if len(i) > 100:       
                    w.write('>ORF_rev_{}\n{}\n'.format(count, ''.join(str(i)))) #write function sooooooooooooo slow don't know why, change to print to see results
                    count+=1
                    
            break
            #REGEX FOR THE WIN
            #print(string, file=filename)

######################################################################## 
######################### DISTANCE MATRIX TOOL #########################
########################################################################    

def distance_matrix():
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
        w.write(d_matrix)
    
if __name__ == '__main__':
    #print(compute_gc())   
    #print(compute_dinucleo()) 
    #print(compute_aa())
    #print(complementDNA())  
    #print(ORF_finder())
    print(distance_matrix())
    

    
    
    
    
    



