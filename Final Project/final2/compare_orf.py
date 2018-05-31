import sys
import re
from Bio import SeqIO

def compare_orf(input_files):
    glimmer = open(input_files[0], 'r')
    data1 = glimmer.readlines()

    orf = open(input_files[1], 'r')
    data2 = orf.readlines()
    
    #proteome: uniprot_03.tab TRANSLATED_ORF
    #genome: 03.fa.txt.longorf ORF\ finder\ FASTA\:\ 03.fa.txt
    correct_glim = [rec.seq for rec in SeqIO.parse(input_files[0],'fasta')]
    correct_orf = [rec.seq for rec in SeqIO.parse(input_files[1],'fasta')]
    #Used Biopython to parse glimmer files
   
    TP=0 # true positive
    TP_list = []
    for line in correct_orf: 
        if line in correct_glim:
            TP+=1 #counts up each TP 
            TP_list.append(line)
    #return TP 
    
    FP=0 # false positive
    FP_list = []
    for line in correct_orf: 
        if  line not in correct_glim:
            FP+=1 #counts up each FP 
            FP_list.append(line)
    #return FP 
    
    FN=0 # false negative
    FN_list = []
    for line in correct_glim: 
        if line not in correct_orf:
            FN+=1 #counts up each FN
            FN_list.append(line)
    return TP,FP,FN 
    
    
def f1_score(): #F1 SCORE
    TP = compare_orf(input_files)[0]
    FP = compare_orf(input_files)[1]
    FN = compare_orf(input_files)[2]
    if TP > 0:
        precision = float(TP)/(TP+FP) # calculates precision with TP/(TP+FP)
        recall = float(TP)/(TP+FN)# calculates recall with TP/(TP+FN)
        return 2*((precision*recall)/(precision+recall)) # f1 score
    else:
        return 0
            
if __name__ == '__main__':
    input_files = sys.argv[1:]
    #print(compare_orf(input_files))
    print(f1_score())

        
    
        
    
