read -p "Which database would you like to extract? (nt/nr):" database

while ! [[ $database == "nt" ]] && ! [[ $database == "nr" ]];
	do
	read -p "Please name the appropriate database for extraction (nt/nr):" database
	wait
done

if [[ $database == "nt" ]]; then
	for FILE in nt.*.tar.gz;
	do
		mkdir ${FILE}_dir
		tar -xvzf ${FILE} -C ./${FILE}_dir/
	done
	
elif [[ $database == "nr" ]]; then
	for FILE in nr.*.tar.gz;
	do
		mkdir ${FILE}_dir
		tar -xvzf ${FILE} -C ./${FILE}_dir/
	done
fi
