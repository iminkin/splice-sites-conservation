cat $1 | while read line
do
	l=($line)
	if [ "${l[0]}.fa.gz" == "$2" ]; then
		tmp_file="$4/${l[1]}.zip"
		datasets download genome accession "${l[1]}" --include=genome --filename="$tmp_file"
		tmp_dir="$4/${l[0]}"
		unzip $tmp_file -d "$tmp_dir"
		data_dir="$tmp_dir/ncbi_dataset/data/${l[1]}"
		fna_file=`ls $data_dir`
		gzip -v3 "$data_dir/$fna_file"
		fna_file=`ls $data_dir`
		mv "$data_dir/$fna_file" "$3/${l[0]}.fa.gz"
		rm $tmp_file
		rm -rf $tmp_dir
	fi
done
