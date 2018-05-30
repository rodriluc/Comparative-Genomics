import sys

def compare_orf(input_files):
    glimmer = open(input_files[0], 'r')
    data1 = glimmer.readlines()
    glimmer_list = data1[1::2]
    
    orf = open(input_files[1], 'r')
    data2 = orf.readlines()
    orf_list = data2[1::2]
    
    for line in glimmer_list:
        if not line.startswith('>'):
            line.strip('\n')
    print(glimmer_list)
    TP=0
    TP_list = []
    for i in orf_list: #TP
        if i in glimmer_list:
            TP+=1
            TP_list.append(i)
    return TP_list
            
if __name__ == '__main__':
    input_files = sys.argv[1:]
    print(compare_orf(input_files))

        
    
        
    
