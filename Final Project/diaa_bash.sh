
for file in *.fasta
do 
echo $file
python2 -c "import diaa2; diaa.compute_diaa('$file')"
done


