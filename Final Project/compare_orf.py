import sys
import re
from Bio import SeqIO

def compare_orf(input_files):
    glimmer = open(input_files[0], 'r')
    data1 = glimmer.readlines()
    #print(data1)
    #glimmer_list = data1[1::2]

    orf = open(input_files[1], 'r')
    data2 = orf.readlines()
    #orf_list = data2[1::2]
    
    #proteome: uniprot_03.tab TRANSLATED_ORF
    #genome: 03.fa.txt.longorf ORF\ finder\ FASTA\:\ 03.fa.txt
    correct_glim = [rec.seq for rec in SeqIO.parse(input_files[0],'fasta')]
    correct_orf = [rec.seq for rec in SeqIO.parse(input_files[1],'fasta')]
    #print(correct_glim)
    '''for line in correct_glim:
        print(line)
        break'''
     
    
    '''templist1 = [] 
    for line in data1:
        if line.startswith('>'):
            nextLine = data1[1]
        #nextLine = data1[line+1]
        templist1.append(nextLine)
    correct_glim = templist1#''.join(templist1)
    return correct_glim 
    
    #correct_glim = open()
   
    for counter in range(0,len(data1)):
         if data1[counter].startswith('>'):
            wholeseq=[]
            counter2=1
            while counter2!=0:
                if data1[counter+counter2].startswith('>')==False:
                    wholeseq.append(data1[counter+counter2])
                    counter2=counter2+1
                else:
                    counter=0
    correct_glim = ''.join(wholeseq)                
    return wholese
                    
    
    templist2 = []
    for line in data2:
        if not line.startswith('>'):
            line = line.strip('\n')
            templist2.extend(line)                  
    correct_orf = ''.join(templist2)
    #return correct_orf
    lines=0
    for line in correct_orf:
        lines+=1
    return lines'''
   
    TP=0 # true positive
    TP_list = []
    for line in correct_orf: #do i need line.strip()
        #nextLine = next(correct)
        if line in correct_glim:
            TP+=1
            TP_list.append(line)
    #return TP 
    
    FP=0 # false positive
    FP_list = []
    for line in correct_orf: 
        #nextLine = next(correct)
        if  line not in correct_glim:
            FP+=1
            FP_list.append(line)
    #return FP 
    
    FN=0 # false negative
    FN_list = []
    for line in correct_glim: 
        #nextLine = next(correct)
        if line not in correct_orf:
            FN+=1
            FN_list.append(line)
    return TP,FP,FN 
    
def sn_sp():
#sn = TP/(TP+FN)
#sp = TP/(TP+FP)
    TP = compare_orf(input_files)[0]
    FP = compare_orf(input_files)[1]
    FN = compare_orf(input_files)[2]
#SENSITIVITY
    sn = TP/(TP+FN)
#SPECIFICITY
    sp = TP/(TP+FP)
    print('sn = ',sn,'\n','sp = ',sp) 
    
def f1_score(): #F1 SCORE
    TP = compare_orf(input_files)[0]
    FP = compare_orf(input_files)[1]
    FN = compare_orf(input_files)[2]
    if TP > 0:
        precision = float(TP)/(TP+FP) # TP/(TP+FP)
        recall = float(TP)/(TP+FN)# TP/(TP+FN)
        return 2*((precision*recall)/(precision+recall))
    else:
        return 0
            
if __name__ == '__main__':
    input_files = sys.argv[1:]
    #print(compare_orf(input_files))
    #print(sn_sp())
    print(f1_score())

        
    
        
    
