#!/usr/bin/env bash
for file in *.fa.txt
do 
echo $file 
python2 grand_finale2.py $file
done

'''for file in *.fasta
do 
echo $file
python2 -c "import grand_finale2; grand_finale2.compute_aa('$file')"
done

for file in *.fasta
do 
echo $file
python2 -c "import grand_finale2; grand_finale2.compute_diaa('$file')"
done

python2 distance_matrix.py 03.fa.txt 28.fa.txt 43.fa.txt 48.fa.txt 50.fa.txt > distance
cat distance
tr -s " " " " < distance > distance1
tr -d "[" < distance1 > distance2
tr -d "]" < distance2 > distance3
sed -i 's/^ *//' distance3
echo -e "03.fa.txt 28.fa.txt 43.fa.txt 48.fa.txt 50.fa.txt\n$(cat distance3)" > distance3
cat distance3

belvu -T R -o tree distance3 > newwick_format_file'''



