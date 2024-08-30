cat $1 | while read line
do
	l=($line)
	if [ "${l[0]}.fa.gz" == "$2" ]; then
		wget -O "$3/$2" "https://dnazoo.s3.wasabisys.com/${l[1]}"/"${l[2]}"
	fi
done
