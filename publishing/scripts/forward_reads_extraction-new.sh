# Bash script used to extract the forwards reads from the raw transposition  paired-end sequencing data (Novogene output).

# We need to know the sequences that contain the sequencing primer as the forward read from the unzipped fasta file. 

## To extract all reads that contain the sequencing read as the forward reads

echo "Changing to data directory"

#cd ../processed-data/

name=$(echo "$(pwd)" | sed 's:.*/::') # name of the folder as the name for the processing files 

echo "Unzipping fq.gz files"
for i in *.fq.gz
do
	gunzip -k ${i}
done
# this will usually output two .fq files, one per each end of the paired-end sequencing data

for i in *.fq
do
	echo " Taking all sequences that contain the sequencing primer while keeping the original format"

	grep --no-group-separator -i -F -B 1 -A 2 "tttaccgaccgttaccgaccgttttcatcccta" $i > ${name}_cleaned_forward_reads_${i}

done


echo " Merging both files together"

## Concatenate all files together
cat ${name}_cleaned* > ${name}_merged_cleaned_forward_reads.fq

echo " Compressing the files for analysis "
### Compressing the files for analysis 
gzip -k ${name}_cleaned*
gzip -k ${name}_merged*

