## script meant to run before the Chromosome_refseq2roman_new_unix.sh

for i in $(find . -name "*_new.wig" -print)
do
	rm $i
done

