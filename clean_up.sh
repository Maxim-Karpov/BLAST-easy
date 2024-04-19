read -p "Which database would you like to clean up (delete the compressed db volumes)? (nt/nr):" database

while ! [[ $database == "nt" ]] && ! [[ $database == "nr" ]];
	do
	read -p "Please name the appropriate database to clean up (nt/nr):" database
	wait
done

for FILE in ${database}.*.tar.gz;
do
rm -r ${FILE}
done

