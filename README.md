BLAST-easy is a collection of Bash scripts which allows for easy installation of the local NCBI full nucleotide (nt) database and conduction of multithreaded, memory-aware BLASTn searches with lower RAM requirements. Currently (14/04/24), the full NCBI nt database requires around 450 GB of disk space when unpacked, thus requiring just as much RAM to run via conventional means. BLAST-easy sets up the nt database and runs a BLAST search in such a way that the database is cached based on user-specified thread and RAM requirements, permitting the BLASTing of the entire NCBI nt database with only 1 core and 16 GB of RAM at user's disposal. The utility of BLAST-easy will increase with time due to the exponential growth of the sequence data on the NCBI database.
<br/>
<br/>
**Requirements (14/04/24 - will increase with the future growth of the database)**
- 1 core
- 16 GB RAM
- 900 GB disk space (~500 GB post-installation; SSD or NVMI storage devices are highly recommended)
<br/>
**Database installation instructions**<br/>
To set up the NCBI nt database, the scripts should be executed inside your local working BLAST search directory in the following order:
  1) download_db.sh - downloads the compressed NCBI nt database by parts.
  2) extract_db.sh - extracts each volume of the downloaded NCBI nt database.
  3) add_metadata_files.sh - adds necessary metadata to each extracted volume of the database for multithreading purposes and significantly lower RAM requirements.
<br/>
<br/>
**BLAST search instructions**<br/>
Run the script BLAST_search.sh in your BLAST-easy directory with a given query FASTA file of choice (e.g. Example_BLAST_query.fasta) in the same directory:

```
BLAST_search.sh Example_BLAST_query.fasta
```

Without additional user input, the program will launch a single BLASTn thread, consuming a single CPU core and under 16 GB of RAM at a time. If higher core and RAM availability is specified, the program will calculate the permissive number of BLASTn threads to run automatically. The number of threads is a slight underestimation of the number of BLASTn processes you may feasibly run, the closer your specifications are to 16 GB of RAM, the more accurate the program's estimation is. This means that at quantities of RAM you may want to specify the -m paramater to be of higher value than the quantity of RAM at your disposal if you would like to optimise the speed of the BLAST search.
<br/>
<br/>
Further usage: BLAST_search.sh [options] {query file name (e.g. query.fasta)}

   -n                    number of cores available to use (default=1)

   -m                    GB of RAM available for use (default=16)

   -max_seqs             BLASTn max_seqs parameter (default=100)

   -max_hsps             BLASTn max_hsps parameter (default=1)

   -e_val                Expect value (E) for saving hits (default=1e-5)

   -outfmt               format of the BLAST output (default='7 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore score pident')

   -h, --help            print description of command line arguments

<br/>
**Future work**<br/>
  1) Ability to use nr database.
  2) Single-script installation of the database of choice.
  3) Ability to run other types of BLAST searches.
  4) Improve BLAST process number estimation algorithm.
