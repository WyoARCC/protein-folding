if [[ $2 == *.fasta ]]
then
	newfile=$(echo $2 | sed 's/_numID//')
	cp $2 $newfile
	list=$(grep -n ">" $2 | awk -F ":" '{print $1}' | tr '\n' ' ')
	
	for j in $list
	do
		line=$(head -n $j $2 | tail -1)
		l2=${line#>seq_}
		number=${l2%_}
		sed -i "s/$line/$(head -n $number $1 | tail -1)/" $newfile
		echo $number >> /project/arcc-students/ecrane3/protfold/anotherfile.txt
	done
	printf "$newfile\n"

else
	newfile=$(echo $2 | sed 's/_num//;s/overview/cazyme_out/;s/prediction_results/signalp_out/')

	cp $2 $newfile
	max=$(sed -n "$=" $newfile)

	for ((j=1;j<=${max};j++))
	do
        	line=$(head -n $j $2 | tail -1)
	        if [[ $line == \seq_* ]]
        	then
                	seq=${line%%	*}
        	        seq=${seq%% *}
	                num=${seq#seq_}
                	sed -i "s/$seq	/$(head -n $num $1 | tail -1)	/" $newfile
        	        sed -i "s/$seq  */$(head -n $num $1 | tail -1)	/" $newfile
	                echo $seq >> /project/arcc-students/ecrane3/protfold/anotherfile.txt
                	echo $num >> /project/arcc-students/ecrane3/protfold/anotherfile.txt
        	fi
	done
	printf "$newfile\n"
fi
