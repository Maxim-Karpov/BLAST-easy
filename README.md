# Overview
BLAST-easy is a collection of Bash scripts which allows for easy installation of the local NCBI full nucleotide (nt) database and conduction of multithreaded, memory-aware BLASTn searches with lower RAM requirements, ideal for personal computer systems, removing the need for Cloud and HPC environments. Currently (14/04/24), the full NCBI nt database takes up around 450 GB of disk space when unpacked, thus requiring just as much RAM to run via conventional means. BLAST-easy sets up the nt database and runs a BLAST search in a manner which allows the database to be cached based on the user-specified thread and RAM requirements, permitting the BLASTing of the entire NCBI nt database with only 1 core and 16 GB of RAM at user's disposal. The utility of BLAST-easy will only become more apparent with time due to the exponential growth of the sequence data on the NCBI database.

## Requirements
- 1 core
- 16 GB RAM
- 900 GB disk space (~500 GB post-installation, will increase with the future growth of the database; SSD or NVMe storage devices are highly recommended)
- BLAST+ installation and the availability of blastn binary in PATH environmental variable

## Database installation instructions
To set up the NCBI nt database, the scripts should be executed inside your local working BLAST-easy directory:
  1) git clone https://github.com/Maxim-Karpov/BLAST-easy
  2) cd ./BLAST-easy
  3) download_db.sh - downloads the compressed NCBI nt database by parts.
  4) extract_db.sh - extracts each volume of the downloaded NCBI nt database.
  5) add_metadata_files.sh - adds necessary metadata to each extracted volume of the database for multithreading purposes and significantly lower RAM requirements.
  6) clean_up.sh (Optional) -  removes the compressed database volumes.

## BLAST search instructions
Run the script BLAST_search.sh in your BLAST-easy directory with a given query FASTA file of choice (e.g. Example_BLAST_query.fasta) in the same directory:

```
bash BLAST_search.sh Example_BLAST_query.fasta
```

Without additional user input, the program will launch a single BLASTn thread, consuming a single CPU core and under 16 GB of RAM at a time. If higher core and RAM availability is specified, the program will automatically calculate the permissive number of BLASTn threads to run. The number of threads is a slight underestimation of the number of BLASTn processes you may feasibly run, the closer your specifications are to 16 GB of RAM, the more accurate the program's estimation is. This means that, when large amounts of RAM are available, you may want to specify the -m paramater to be of a higher value than the actual quantity of RAM at your disposal if you would like to optimise the speed of the BLAST search (i.e. run more simultaneous BLASTn processes). The automatic process number assignment algorithm however will make sure that "out of memory" errors are not encountered.

```
Further usage: BLAST_search.sh {query file name (e.g. query.fasta)} [options] 

   -n                    number of cores available for use (default=1)
   -m                    GB of RAM available for use (default=16)
   -max_seqs             BLASTn max_seqs parameter (default=100)
   -max_hsps             BLASTn max_hsps parameter (default=1)
   -e_val                Expect value (E) for saving hits (default=1e-5)
   -outfmt               format of the BLAST output (default='7 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore score pident')
   -h, --help            print description of command line arguments
```

## Future implementations
  1) Ability to use the full NCBI nr protein sequence database.
  2) Single-script installation of the database of choice.
  3) Ability to run other types of BLAST searches.
  4) Improved BLAST process number estimation algorithm.
  5) Smarter BLAST output merging from different database volumes.
