for FILE in nt.*.tar.gz;
do
mkdir ${FILE}_dir
tar -xvzf ${FILE} -C ./${FILE}_dir/
done
