#!/bin/bash
counter=0
for FILE in nt.*.tar.gz;
do
cd ${FILE}_dir
rm nt.nal
touch nt.nal

if [[ "$counter" -le 9 ]]; then
printf "#\n# Alias file created: Apr 1, 2024  7:29 PM\n#\nTITLE Nucleotide collection (nt)\nDBLIST nt.00${counter}" >> nt.nal
counter=$((counter+1))	
elif [[ "$counter" -le 99 ]]; then
printf "#\n# Alias file created: Apr 1, 2024  7:29 PM\n#\nTITLE Nucleotide collection (nt)\nDBLIST nt.0${counter}" >> nt.nal
counter=$((counter+1))
else
printf "#\n# Alias file created: Apr 1, 2024  7:29 PM\n#\nTITLE Nucleotide collection (nt)\nDBLIST nt.${counter}" >> nt.nal
counter=$((counter+1))
fi
cd ..
ln -s $(pwd)/nt.000.tar.gz_dir/nt.ndb $(pwd)/${FILE}_dir/nt.ndb
ln -s $(pwd)/nt.000.tar.gz_dir/nt.nos $(pwd)/${FILE}_dir/nt.nos
ln -s $(pwd)/nt.000.tar.gz_dir/nt.ntf $(pwd)/${FILE}_dir/nt.ntf
ln -s $(pwd)/nt.000.tar.gz_dir/nt.not $(pwd)/${FILE}_dir/nt.not
ln -s $(pwd)/nt.000.tar.gz_dir/nt.nto $(pwd)/${FILE}_dir/nt.nto
ln -s $(pwd)/nt.000.tar.gz_dir/taxdb.btd $(pwd)/${FILE}_dir/taxdb.btd
ln -s $(pwd)/nt.000.tar.gz_dir/taxdb.bti $(pwd)/${FILE}_dir/taxdb.bti
ln -s $(pwd)/nt.000.tar.gz_dir/taxonomy4blast.sqlite3 $(pwd)/${FILE}_dir/taxonomy4blast.sqlite3
done 
