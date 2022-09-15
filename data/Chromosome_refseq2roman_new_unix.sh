
############ erase first all the "_roman.wig" already created files , when running this script. Otherwise it will create another wig file on top of that##########################

#for i in $(find . -name "*_roman.wig" -print)
#do
#rm $i
#done


# for loop over the wig files in all subdirectories 

for i in $(find . -name "*.wig" -print)  #  list all wig files in the subdirectories of this root folder 
do

	sed 's/variablestep/variableStep/g' ${i} > ${i}_new.wig

done



