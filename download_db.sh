read -p "Which database would you like to download? (nt/nr):" database

if [[ $database == "nt" ]] || [[ $database == "nr" ]]; then
	echo "Please select the appropriate database (nt/nr):"
	exit 1
fi

if [[ $database == "nt" ]]; then
	wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
	
elif [[ $database == "nr" ]]; then
	wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz"
fi
