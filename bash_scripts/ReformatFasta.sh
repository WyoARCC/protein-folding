cp $1 $2
list=$(grep -n ">" $2 | awk -F ":" '{print $1}' | tr '\n' ' ')

for j in $list
do
	line=$(head -n $j $1 | tail -1)
	sed -i "s/$line/>seq_$j/" $2
done
