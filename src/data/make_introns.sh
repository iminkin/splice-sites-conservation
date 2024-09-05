for transcript in $1/transcripts/$2/*
do
	IFS='/'; p=($transcript)
	chr=${p[5]}; transcript_id=${p[6]}; IFS=' '
	file=$chr/$transcript_id
	if [ ! -d "$1/introns/$chr" ]; then
		mkdir $1/introns/$chr
	fi
	bedtools subtract -a $1/transcripts/$file -b $1/exons/$file > $1/introns/$file
done

cat $1/introns/$2/* > $1/introns_all/$2
