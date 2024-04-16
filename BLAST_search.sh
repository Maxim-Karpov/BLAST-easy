#!/bin/bash

#Usage description
help_present=0
if [[ $* == *"-h"* ]] || [[ $* == *"--help"* ]]; then
	help_present=1
fi
if [[ $help_present == 1 ]]; then
	echo "Usage: $0 {query file name (e.g. query.fasta)} [options] "
	echo
	echo "   -n                    number of cores available to use (default=1)"
	echo "   -m                    GB of RAM available for use (default=16)"
	echo "   -max_seqs             BLASTn max_seqs parameter (default=100)"
	echo "   -max_hsps             BLASTn max_hsps parameter (default=1)"
	echo "   -e_val                Expect value (E) for saving hits (default=1e-5)"
	echo "   -outfmt               format of the BLAST output (default='7 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore score pident')"
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
		echo $3
		OUTFMT=$3
		shift 2
fi
done

#Set variable names/values and defaults
C="${n:-1}"
RAM="${m:-16}"
MAX_SEQS="${max_seqs:-100}"
MAX_HSPS="${max_hsps:-1}"
E_VAL="${e_val:-1e-5}"

if [ $OUTFMT -z ]; 
then
	OUTFMT=("7 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore score pident")
fi

total=0
file_number=0

#Calculate size of database and number of database directories
for FILE in *_dir/;
do
	sz=$(du -hb ${FILE} | cut -f1)
	total=$((total + sz))
	file_number=$((file_number+1))
done

#Calculate process limit assigned by available RAM
sz_zero=$(du -hb nt.000.tar.gz_dir | cut -f1)
average_file_size=$((total/file_number/1000000000))
RAM_process_limit=$(($RAM/((($average_file_size+1)+($sz_zero/1000000000))/2)))
RAM_process_limit=$((RAM_process_limit-1))
ar=($RAM_process_limit $C)

#Calculate maximum processes based on available cores and RAM
max_processes=${ar[0]}
for n in "${ar[@]}";
do
	((n < max_processes)) && max_processes=$n
done

echo "Maximum number of processes:" $max_processes

#Run BLAST
counter=0
for FILE in *_dir/;
do
	ln -s $(pwd)/$QUERY_NAME $(pwd)/$FILE/$QUERY_NAME
	cd $FILE
	blastn -query $QUERY_NAME -db nt -max_target_seqs $MAX_SEQS -max_hsps $MAX_HSPS -evalue $E_VAL -outfmt "$OUTFMT" -out out_${counter}.txt -num_threads 1 &
	counter=$((counter+1))
	cd ..
	if ! (($counter % $max_processes)); then
		wait
fi
done

#Kill any running BLAST processes in case of a manual interrupt
kill_blastn() {
	trap SIGINT
	pkill -f blastn
}
trap "kill_blastn" SIGINT 

#Compile all database thread outputs into one
touch out.txt
for FILE in *_dir/;
do
	echo $(pwd)/$FILE/out_*.txt >> out.txt
done
