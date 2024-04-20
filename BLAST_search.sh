#!/bin/bash

#Usage description
help_present=0
if [[ $* == *"-h"* ]] || [[ $* == *"--help"* ]]; then
	help_present=1
fi
if [[ $help_present == 1 ]]; then
	echo "Usage: $0 {query file name (e.g. query.fasta)} [options] "
	echo
	echo "   -n                    number of cores available to use [default=1]"
	echo "   -m                    GB of RAM available for use [default=5]"
	echo "   -db                   database selection (nt/nr, corresponding to blastn/blastp search, respectively) [default=nt]"
	echo "   -max_seqs             BLAST max_seqs parameter [default=100]"
	echo "   -max_hsps             BLAST max_hsps parameter [default=1]"
	echo "   -e_val                Expect value (E) for saving hits [default=1e-5]"
	echo "   -outfmt               format of the BLAST output [default='7 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore score pident']"
	echo "   -rc                   continuously refresh cache during BLAST search (requires user privileges) [0/1] [default=0]"
	echo "   -custom               custom BLAST options (e.g. '-perc_identity 20 -gapopen 10 -gapextend 5 -html') [default=]"
	echo "   -h, --help            print description of command line arguments"
	echo
	exit 1
fi

#Check for mandatory argument
if [[ -z "$1" ]]; then
	echo "BLAST query file name must be specified (e.g. query.fasta)."
	exit 1
fi

#Process arguments and options
QUERY_NAME=$1
while ! [ -z "$2" ]
do
	if [ "$2" = "-n" ];
	then
		n=$3
		shift 2
	elif [ "$2" = "-m" ];
	then
		m=$3
		shift 2
	elif [ "$2" = "-db" ];
	then
		db=$3
		shift 2
	elif [ "$2" = "-max_seqs" ];
	then
		max_seqs=$3
		shift 2
	elif [ "$2" = "-max_hsps" ];
	then
		max_hsps=$3
		shift 2
	elif [ "$2" = "-e_val" ];
	then
		e_val=$3
		shift 2
	elif [ "$2" = "-outfmt" ];
	then
		OUTFMT=$3
		shift 2
	elif [ "$2" = "-rc" ];
	then
		rc=$3
		shift 2
	elif [ "$2" = "-custom" ];
	then
		CUSTOM=$3
		shift 2
fi
done

#Set variable names/values and defaults
C="${n:-1}"
RAM="${m:-5}"
DB="${db:-nt}"
MAX_SEQS="${max_seqs:-10}"
MAX_HSPS="${max_hsps:-1}"
E_VAL="${e_val:-1e-5}"
RC="${rc:-0}"

if [ $OUTFMT -z ]; 
then
	OUTFMT=("7 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore score pident")
fi

if [ $CUSTOM -z ]; 
then
	CUSTOM=
fi

total=0
file_number=0

#Calculate size of the first database volume and number of database directories
#for FILE in *_dir/;
#do
#	sz=$(du -hb ${FILE} | cut -f1)
#	total=$((total + sz))
#	file_number=$((file_number+1))
#done

#Calculate size of nsq file and number of database directories
for FILE in *_dir/;
do
	cd $FILE
	for inFILE in *.nsq; do
		sz=$(du -hb ${inFILE} | cut -f1)
	done
	if [[ $sz -ge $total ]]; then
		total=$sz
	fi
	file_number=$((file_number+1))
	cd ..
done

if [[ $file_number -ge 101 ]]; then
numbers="000"
else
numbers="00"
fi

#Calculate process limit assigned by available RAM, based on average file size and the size of the first volume
#sz_zero=$(du -hb ${DB}.${numbers}.tar.gz_dir | cut -f1)
#average_file_size=$((total/file_number/1000000000))
#RAM_process_limit=$(($RAM/((($average_file_size+1)+($sz_zero/1000000000))/2)))
#RAM_process_limit=$((RAM_process_limit-1))
#ar=($RAM_process_limit $C)


#Calculate process limit assigned by available RAM, based on the maximal size of the .nsq file
total1=$(($total / 1000000000 + 2))
RAM_process_limit=$(($RAM/$total1))
ar=($RAM_process_limit $C)


#Calculate maximum processes based on available cores and RAM
max_processes=${ar[0]}
for n in "${ar[@]}";
do
	((n < max_processes)) && max_processes=$n
done

echo "Maximum number of processes:" $max_processes

if [[ $max_processes -le 0 ]]; then
	echo "Because the calculated maximum number of processes is lower than 1, the program may function at a slower rate, or may not function at all (e.g. crash due to an 'Out Of Memory' error), which will be evident during the caching of the first database volume"
	answer_t=
	while [[ $answer_t != "Y" ]]; do
		read -p "Would you still like to proceed? (Y/N)" answer_t
		wait
		if [[ $answer_t == "N" ]]; then
			exit 1
		fi
	done
	max_processes=1
fi

#Run BLAST
counter=0
if [[ $DB == "nt" ]]; then
	for FILE in *_dir/;
	do
		if [[ $RC == 1 ]]; then
		sync && echo 3 | sudo tee /proc/sys/vm/drop_caches>/dev/null 
		fi
		
		if [[ -f $(pwd)/$FILE/$QUERY_NAME ]]; then
			rm $(pwd)/$FILE/$QUERY_NAME
			ln -s $(pwd)/$QUERY_NAME $(pwd)/$FILE/$QUERY_NAME
		else
			ln -s $(pwd)/$QUERY_NAME $(pwd)/$FILE/$QUERY_NAME
		fi
		
		cd $FILE
		blastn -query $QUERY_NAME -db $DB -max_target_seqs $MAX_SEQS -max_hsps $MAX_HSPS -evalue $E_VAL -outfmt "$OUTFMT" -out out_${counter}.txt -num_threads 1 $CUSTOM &
		counter=$((counter+1))
		cd ..
		if ! (($counter % $max_processes)); then
			wait
		fi
	done
	
elif [[ $DB == "nr" ]]; then
	for FILE in *_dir/;
	do

		if [[ $RC == 1 ]]; then
		sync && echo 3 | sudo tee /proc/sys/vm/drop_caches>/dev/null 
		fi
		
		if [[ -f $(pwd)/$FILE/$QUERY_NAME ]]; then
			rm $(pwd)/$FILE/$QUERY_NAME
			ln -s $(pwd)/$QUERY_NAME $(pwd)/$FILE/$QUERY_NAME
		else
			ln -s $(pwd)/$QUERY_NAME $(pwd)/$FILE/$QUERY_NAME
		fi
			
		cd $FILE
		blastp -query $QUERY_NAME -db $DB -max_target_seqs $MAX_SEQS -max_hsps $MAX_HSPS -evalue $E_VAL -outfmt "$OUTFMT" -out out_${counter}.txt -num_threads 1 $CUSTOM &
		counter=$((counter+1))
		cd ..
		if ! (($counter % $max_processes)); then
			wait
		fi
	done
fi

#Kill any running BLAST processes in case of a manual interrupt
kill_blastn() {
	trap SIGINT
	pkill -f blastn
}
trap "kill_blastn" SIGINT 

#Compile all database thread outputs into one
rm out.txt
touch out.txt
for FILE in *_dir/;
do
	cat $(pwd)/$FILE/out_*.txt >> out.txt
done

#filter no match entries and parse files for ease of analysis
touch out.filtered.txt
touch out.nomatch.txt
touch out.hits_only.rm1.txt


MAX_SEQS=10
OUTFMT=("7 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore score pident")

counter=1
while [[ $counter -le $MAX_SEQS ]];
do
cat ./out.txt | grep "${counter} hits" -A ${counter} -B 3 >> out.filtered.txt 
counter=$((counter+1))
done

cat ./out.txt | grep "0 hits" -B 2 >> out.nomatch.txt

counter=1
while [[ $counter -le $MAX_SEQS ]];
do
cat ./out.filtered.txt | grep "${counter} hits" -A ${counter} >> out.hits_only.rm1.txt
counter=$((counter+1))
done
cat out.hits_only.rm1.txt | grep -v " hits found" > out.hits_only.rm2.txt
cat out.hits_only.rm2.txt | grep -v -e "--" > out.hits_only.txt
rm out.hits_only.rm1.txt
rm out.hits_only.rm2.txt

#add a header to the BLAST result
phrase=
for WORD in $OUTFMT;
do
phrase=${phrase}${WORD}"\t"
done
OUTPUT="$(echo $OUTFMT | awk '$1=$1' OFS="\t")"
header="$(echo "$OUTPUT" | cut -f 2- -d '	')"
echo $header | cat - out.hits_only.txt > temp && mv temp out.hits_only.txt
tr " " "\t" < out.hits_only.txt > temp && mv temp out.hits_only.txt
