cat $1 | while read line
do
	l=($line)
	if [ "${l[0]}.fa.gz" == "$2" ]; then
		tmp_file="$4/${l[1]}.zip"
#		datasets download genome accession "${l[1]}" --include=genome --filename="$tmp_file"
#		unzip $tmp_file -d "$4"
		data_dir="$4/ncbi_dataset/data/${l[1]}"
		fna_file=`ls $data_dir`
		gzip "$data_dir/$fna_file"
	fi
done
