#!/bin/bash

#Usage description
help_present=0

if [[ $* == *"-h"* ]] || [[ $* == *"--help"* ]]; then
	help_present=1
fi

if [[ $help_present == 1 ]]; then
	echo "Usage: $0 [options] "
	echo
	echo "   -db                   database selection [nt/nr]"
	echo "   -checkpoint           start the installation from a certain process [1/2/3] [default=1]: 
	                  1) Download database 
	                  2) Decompress volumes 
	                  3) Add metadata files "
	echo "   -rc                   continuously refresh cache during the installation process search (requires user privileges) [0/1] [default=0]"
	echo "   -h, --help            print description of command line arguments"
	echo
	exit 1
fi


#Process arguments and options
while ! [ -z "$1" ]
do
	if [ "$1" = "-db" ];
	then
		db=$2
		shift 2
	elif [ "$1" = "-checkpoint" ];
	then
		cp=$2
		shift 2
	elif [ "$1" = "-rc" ];
	then
		rc=$2
		shift 2
fi
done


#User inputs
if [[ $db != "nr" ]] && [[ $db != "nt" ]]; then
	read -p "Which database would you like to install? (nt/nr):" database

	while ! [[ $database == "nt" ]] && ! [[ $database == "nr" ]];
		do
		read -p "Please name the appropriate database to install (nt/nr):" database
		wait
	done
fi

#Set variable names/values and defaults
C="${cp:-1}"
RC="${rc:-0}"


#Download compressed database
if [[ $C == 1 ]] || [[ $C == "" ]]; then
	echo "Checkpoint 1: Downloading ${database} database"
	wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/${database}.*.tar.gz"
	C=2
fi


#Unpack compressed database files
if [[ $C == 2 ]] || [[ $C == "" ]]; then
	echo "Checkpoint 2: Extracting ${database} database volumes"
	counter=0
	for FILE in ${database}.*.tar.gz;
	do
	
	#Cache refresh option
	if [[ $RC == 1 ]]; then
	sync && echo 3 | sudo tee /proc/sys/vm/drop_caches>/dev/null 
	fi

	counter=$((counter+1))

	mkdir ${FILE}_dir
	tar -xvzf ${FILE} -C ./${FILE}_dir/
	done
	C=3
fi


#Enrich with metadata
if [[ $C == 3 ]] || [[ $C == "" ]]; then
	echo "Checkpoint 3: Enriching ${database} database volumes with metadata"
	if [[ $database == "nt" ]]; then

	counter=0
	for FILE in ${database}.*.tar.gz;
	do
	cd ${FILE}_dir

	#Cache refresh
	if [[ -f ${database}.nal ]]; then
		rm ${database}.nal
	fi

	touch ${database}.nal

	if [[ "$counter" -le 9 ]]; then
		printf "#\n# Alias file created: Apr 1, 2024  7:29 PM\n#\nTITLE Nucleotide collection (${database})\nDBLIST ${database}.00${counter}" >> ${database}.nal
		counter=$((counter+1))	
	elif [[ "$counter" -le 99 ]]; then
		printf "#\n# Alias file created: Apr 1, 2024  7:29 PM\n#\nTITLE Nucleotide collection (${database})\nDBLIST ${database}.0${counter}" >> ${database}.nal
		counter=$((counter+1))
	else
		printf "#\n# Alias file created: Apr 1, 2024  7:29 PM\n#\nTITLE Nucleotide collection (${database})\nDBLIST ${database}.${counter}" >> ${database}.nal
		counter=$((counter+1))
	fi
	cd ..
	ln -s -f $(pwd)/${database}.000.tar.gz_dir/${database}.ndb $(pwd)/${FILE}_dir/${database}.ndb
	ln -s -f $(pwd)/${database}.000.tar.gz_dir/${database}.nos $(pwd)/${FILE}_dir/${database}.nos
	ln -s -f $(pwd)/${database}.000.tar.gz_dir/${database}.ntf $(pwd)/${FILE}_dir/${database}.ntf
	ln -s -f $(pwd)/${database}.000.tar.gz_dir/${database}.not $(pwd)/${FILE}_dir/${database}.not
	ln -s -f $(pwd)/${database}.000.tar.gz_dir/${database}.nto $(pwd)/${FILE}_dir/${database}.nto
	ln -s -f $(pwd)/${database}.000.tar.gz_dir/taxdb.btd $(pwd)/${FILE}_dir/taxdb.btd
	ln -s -f $(pwd)/${database}.000.tar.gz_dir/taxdb.bti $(pwd)/${FILE}_dir/taxdb.bti
	ln -s -f $(pwd)/${database}.000.tar.gz_dir/taxonomy4blast.sqlite3 $(pwd)/${FILE}_dir/taxonomy4blast.sqlite3
	done

	elif [[ $database == "nr" ]]; then

	counter=0
	for FILE in ${database}.*.tar.gz;
	do
	cd ${FILE}_dir

	if [[ -f ${database}.pal ]]; then
		rm ${database}.pal
	fi

	touch ${database}.pal

	if [[ "$counter" -le 9 ]]; then
		printf "#\n# Alias file created 04/11/2024 04:54:09\n#\nTITLE All non-redundant GenBank CDS translations+PDB+SwissProt+PIR+PRF excluding environmental samples from WGS projects\nDBLIST ${database}.0${counter}\nNSEQ 721436858\nLENGTH 278554112155" >> ${database}.pal
		counter=$((counter+1))	
	elif [[ "$counter" -le 99 ]]; then
		printf "#\n# Alias file created 04/11/2024 04:54:09\n#\nTITLE All non-redundant GenBank CDS translations+PDB+SwissProt+PIR+PRF excluding environmental samples from WGS projects\nDBLIST ${database}.${counter}\nNSEQ 721436858\nLENGTH 278554112155" >> ${database}.pal
		counter=$((counter+1))
	else
		printf "#\n# Alias file created 04/11/2024 04:54:09\n#\nTITLE All non-redundant GenBank CDS translations+PDB+SwissProt+PIR+PRF excluding environmental samples from WGS projects\nDBLIST ${database}.${counter}\nNSEQ 721436858\nLENGTH 278554112155" >> ${database}.pal
		counter=$((counter+1))
	fi
	cd ..
	ln -s -f $(pwd)/${database}.00.tar.gz_dir/${database}.pdb $(pwd)/${FILE}_dir/${database}.pdb
	ln -s -f $(pwd)/${database}.00.tar.gz_dir/${database}.pos $(pwd)/${FILE}_dir/${database}.pos
	ln -s -f $(pwd)/${database}.00.tar.gz_dir/${database}.ptf $(pwd)/${FILE}_dir/${database}.ptf
	ln -s -f $(pwd)/${database}.00.tar.gz_dir/${database}.pot $(pwd)/${FILE}_dir/${database}.pot
	ln -s -f $(pwd)/${database}.00.tar.gz_dir/${database}.pto $(pwd)/${FILE}_dir/${database}.pto
	ln -s -f $(pwd)/${database}.00.tar.gz_dir/taxdb.btd $(pwd)/${FILE}_dir/taxdb.btd
	ln -s -f $(pwd)/${database}.00.tar.gz_dir/taxdb.bti $(pwd)/${FILE}_dir/taxdb.bti
	ln -s -f $(pwd)/${database}.00.tar.gz_dir/taxonomy4blast.sqlite3 $(pwd)/${FILE}_dir/taxonomy4blast.sqlite3
	done
	fi
fi

