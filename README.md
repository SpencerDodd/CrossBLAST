# CrossBLAST

This is a rough script that I have been using to do genomic research with automated BLAST queries. If you've stumbled on this page, please don't judge the code...I know it's gross and some of it is pretty spaghetti. I'm working on cleaning and fixing it up now, but that's taking a backseat to data generation at the moment.

## Overview

Here's a quick rundown of what each script does, and future functionality that it will have in coming iterations.

### blast_accession.py

This is the workhouse of the scripts. And the most spaghetti. It is run with the following command and arguments:

--------------------------------------------------------
```
python blast_accession.py INPUT_TYPE DATABASE QUERY_SPECIES_NAME QUERY_SUBSPECIES_NAME QUERY_ACCESSION_NUMBER
```

where...

INPUT_TYPE is the style in which the query data will be given to the script
	a: accession number
	f: FASTA file
	cross: cross_blast request (not applicable for shell use, but will be elaborated later)

DATABASE is any of the NCBI's publicly available databases 
	human_genomic
	nt
	refseq_genomic
	etc...

QUERY_SPECIES is the name of the species of the organism the query sequence belongs to
	This is used due to the fact that the NCBI's API only returns phylogenetic information up to the genus level. This value is used to regex match with hits to determine level of relatedness past the genus level

QUERY_SUBSPECIES is the name of the subspecies of the organism the query sequence belongs to
	This is used for the same reason as QUERY_SPECIES

QUERY_ACCESSION_NUMBER is the accession number of the query sequence
	This is only a necessary input if you are using an accession number based query. (and in cross_blast queries)
--------------------------------------------------------

This script queries the NCBI's BLAST servers to look for sequence hits in the inputted database using the blastn algorithm. It outputs the hit information in .csv format delineated by the Family, SubFamily, Genus, Species, and Subspecies phylogenetic levels. 

One important facet of the data analysis is that the % difference field in the output is actually the % difference between the query and hit / 2. This finds the distance to closest common ancestor between the hit and query.
