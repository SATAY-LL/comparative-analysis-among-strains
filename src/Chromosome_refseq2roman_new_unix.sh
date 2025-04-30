
############ erase first all the "_roman.wig" already created files , when running this script. Otherwise it will create another wig file on top of that##########################

#for i in $(find . -name "*_roman.wig" -print)
#do
#rm $i
#done


# for loop over the wig files in all subdirectories 

for i in $(find . -name "*.wig" -print)  #  list all wig files in the subdirectories of this root folder 
do

	sed  '
	s/ref|NC_001133|/I/ 
	s/ref|NC_001134|/II/
	s/ref|NC_001135|/III/  
	s/ref|NC_001136|/IV/
	s/ref|NC_001137|/V/  
	s/ref|NC_001138|/VI/
	s/ref|NC_001139|/VII/  
	s/ref|NC_001140|/VIII/
	s/ref|NC_001141|/IX/  
	s/ref|NC_001142|/X/
	s/ref|NC_001143|/XI/
	s/ref|NC_001144|/XII/  
	s/ref|NC_001145|/XIII/
	s/ref|NC_001146|/XIV/
	s/ref|NC_001147|/XV/  
	s/ref|NC_001148|/XVI/
	s/VariableStep/variableStep/
	' ${i} > ${i}_roman.wig

done



