__author__ = 'Spencer Dodd' # I almost don't want to take credit for this pasta

# NOTE: It's spaghetti. God help your soul if you need to debug or maintain this.

import Tkinter, tkFileDialog
import datetime
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SearchIO
from Bio import Entrez
import xlwt
import xlrd
import csv
import xml.etree.ElementTree as ET
import sys


# -------------------------------------------
# ------ ENSURE IS CORRECT BEFORE RUN! ------
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# ------------ GLOBAL VARIABLES -------------
# -------------------------------------------
# email parameter for NCBI issues with access
Entrez.email = 'dodd.s@husky.neu.edu'
#Entrez.email = 'spencer.evan.dodd@gmail.com'

queries = []
xml_results = []

# holds important information for query hits
	# 0: Accession number
	# 1: % identity
	# 2: phylogeny
	# 3: FASTA sequence (full sequence)
hit_accessions = []
hit_seqs = []
hit_full_seqs = []
hit_percent_identities = []

# date and time information for writing to file
today = datetime.datetime.today().strftime('%y_%m_%d')
hour = datetime.datetime.today().time().hour
minute = datetime.datetime.today().time().minute
second = datetime.datetime.today().time().second

# PHYLOGENETIC DEFINITIONS
# Because the XML returns from qblast are info-sparse, we need to look at 
# phylogenetic info and comparisons in weird ways

Families = [
	'Anguillidae', 			# eel
	'Cervidae',				# deer
	'Bathyergidae', 		# naked mole rat
	'Ursidae',				# asian black bear
	'Salmonidae',			# salmon
	'Phasianidae',			# red junglefowl
	'Hominidae',			# human / hominids
	'Muridae',				# Mouse
	'Balaena'				# Bowhead whale
]
Subfamilies = [
	'Cervinae',				# deer
	'Odocoileinae',			# deer
	'Muntiacinae',			# deer
	'Capreolinae',			# deer
	'Procervulinae',		# deer
	'Heterocephalinae', 	# NMR
	'Salmoninae',			# salmon
	'Phasianidae',			# red junglefowl
	'Murinae'				# Mouse
]
# -------------------------------------------

# allows user to select files to be processed for blast hits
def select_files():

	print '\nSelect FASTA files to BLAST ...'

	# query by fasta file (opened with tkFileDialog)
	root = Tkinter.Tk()
	root.withdraw()
	file_paths = tkFileDialog.askopenfilenames()

	for path in file_paths:

		queries.append(path)

# processes each of the individual queries stored in the queries array from
# 'select_files()'
def handle_requests(query_database, input_type, query_accession):

	print '\nQuerying NCBI database ({0}) ...'.format(query_database)

	# output file to folder "BLAST_QUERIES", 2 levels back from pwd
	if input_type == 'accession':
		
		file_path = '{0}BLAST_QUERIES/'.format(result_output_directory(os.getcwd()))

		# make dir if it doesn't already exist
		if not os.path.exists(file_path):

			os.makedirs(file_path)

	# if input_type is cross
	elif input_type == 'cross':

		file_path = '{0}BLAST_QUERIES/'.format(output_file_path)

		# make dir if it doesn't already exist
		if not os.path.exists(file_path):

			os.makedirs(file_path)

	if input_type == 'fasta':

		attempt_num = 1
		failed_attempt = True
		file_name = '{0}query_(attempt_{1})_({2}h{3}m{4}s).xml'.format(file_path, attempt_num, hour, minute, second)

		while failed_attempt:

			for index, query in enumerate(queries):

				print '\nQuerying for sequence {0}: {1}'.format((index + 1), query)

				print '\nAttempt {0}'.format(attempt_num)

				fasta_string = open(query).read()

				fasta_result_handle = NCBIWWW.qblast(program='blastn', database=query_database, sequence=fasta_string)

				# removes any previous XML results
				del xml_results[:]

				# saves the xml result of the BLAST query
				save_xml_file(file_name, fasta_result_handle)

				# failed_attempt is true if the BLAST result returned a CPU usage error
				# a non-failed attempt should break loop
				failed_attempt = cpu_usage_error(file_name)

				attempt_num += 1

		# stores the file name for xml parsing later (only if query was successful)
		xml_results.append(file_name)

	elif input_type == 'accession': # TODOTODOTODOTODO

		attempt_num = 1
		failed_attempt = True
		file_name = '{0}query_(attempt_{1})_({2}h{3}m{4}s).xml'.format(file_path, attempt_num, hour, minute, second)

		while failed_attempt:

			handle = Entrez.efetch(db='nucleotide', id=query_accession, rettype='fasta')

			record = handle.read()

			print '\nQuerying for sequence {0}: '.format(record.splitlines()[0])

			print '\nAttempt {0}'.format(attempt_num)

			fasta_string = record

			fasta_result_handle = NCBIWWW.qblast(program='blastn', database=query_database, sequence=fasta_string)

			# removes any previous XML results
			del xml_results[:]

			# saves the xml result of the BLAST query
			save_xml_file(file_name, fasta_result_handle)

			# failed_attempt is true if the BLAST result returned a CPU usage error
			# a non-failed attempt should break loop
			failed_attempt = cpu_usage_error(file_name)

			attempt_num += 1

		# stores the file name for xml parsing later (only if query was successful)
		xml_results.append(file_name)

	elif input_type == 'cross':

		attempt_num = 1
		failed_attempt = True
		file_name = '{0}query_(attempt_{1})_({2}h{3}m{4}s).xml'.format(file_path, attempt_num, hour, minute, second)

		while failed_attempt:

			handle = Entrez.efetch(db='nucleotide', id=query_accession, rettype='fasta')

			record = handle.read()

			print '\nQuerying for sequence {0}: '.format(record.splitlines()[0])

			print '\nAttempt {0}'.format(attempt_num)

			fasta_string = record

			fasta_result_handle = NCBIWWW.qblast(program='blastn', database=query_database, sequence=fasta_string)

			# removes any previous XML results
			del xml_results[:]

			# saves the xml result of the BLAST query
			save_xml_file(file_name, fasta_result_handle)

			# failed_attempt is true if the BLAST result returned a CPU usage error
			# a non-failed attempt should break loop
			failed_attempt = cpu_usage_error(file_name)

			attempt_num += 1

		# stores the file name for xml parsing later (only if query was successful)
		xml_results.append(file_name)

# checks to see if global variable hit_accessions is not filled with a CPU usage limit error
# returns true if the XML file it is fed contains the error
def cpu_usage_error(file_name):

	this_hour = datetime.datetime.today().time().hour
	this_minute = datetime.datetime.today().time().minute
	this_second = datetime.datetime.today().time().second

	tree = ET.parse(file_name)
	root = tree.getroot()

	# parses the query's database number from the xml file
	query_db = root[8][0][5][0][0].text

	# if there is a CPU error
	if len(root[8][0]) == 7:

		print '\n---- BLAST CPU ERROR ({0}h {1}m {2}s) ----'.format(this_hour, this_minute, this_second)

		print '\nRe-trying BLAST query ...'

		return True

	# if there was no database reply from the algorithm, but query wasn't rejected
	# (Empty results)
	# Do not retry
	elif len(root[8][0]) == 6 and query_db == '0':

		print '\n---- Other BLAST ERROR (No Results for query) ({0}h {1}m {2}s) ----'.format(this_hour, this_minute, this_second)

		return False

	# if the query went successfully
	else:

		return False

# saves the current XML query object to a file
def save_xml_file(file_name, query_object):

	print '\nSaving XML results ...'

	save_file = open(file_name, 'w')
	save_file.write(query_object.read())
	save_file.close()
	query_object.close()

# parses saved XML data into individual hit data
	# FASTA sequence
	# Phylogeny
	# % identity to query seq
def parse_queries(query_database, query_species, query_subspecies, query_type):

	print '\nparsing ...'

	for file in xml_results:

		print '\nopening file ...'
		result_handle = open(file)

		print '\nXML processing ...'
		blast_record = NCBIXML.read(result_handle)

		query_name = ''

		hit_num = 0

		blast_qresult = SearchIO.read(file, 'blast-xml')
		for hit in blast_qresult:

			for hsp in hit:

				if hit_num == 0:

						query_name = hit.description

				# calculate the % identity of the hit to the query sequence
				hit_identity = ((1.0 * hsp.ident_num / hsp.aln_span) * 100)

				# add accession numbers to the hit_accessions array
				hit_accessions.append( [hit.accession, hit_identity] )

				hit_num += 1
		
		# remove duplicate hits from hit_accessions
		print '\nConsolidating duplicate hits ...'
		
		remove_duplicate_hits(hit_accessions)

		print '\nDone'

		# for each accession number in the hit_accessions:
			# get the FASTA sequence
			# get the phylogeny

		print '\nQuerying GenBank ...'

		all_accessions, full_hit_seqs = query_genbank(query_name, query_database)

		print '\nOutputting Results ...'

		output_full_hits(query_name, full_hit_seqs, all_accessions)

		# remove hits that are not full sequences (sequence length < 16000)
		print '\nRemoving incomplete sequences ...'

		# ----------------------------------------- SORTS HITS BY PHYLOGENY
		# -----------------------------------------------------------------
		# sorts the hits by their phylogenetic relation to query
		# outputs final hit and query data to file directory
		print '\nSorting hits by Phylogeny ...'
		sort_by_phylogeny(query_species, query_subspecies, query_database, query_type)
		# -----------------------------------------------------------------

# sorts the unique hits by phylogeny and outputs the results in
		# Family, Subfamily, Genus, Species determinations of relationship
		# with filenames beginning with % divergence from the query sequence
		#		divergence = (100 - identity) / 2
def sort_by_phylogeny(query_species, query_subspecies, query_database, query_type):

	Family = []
	Subfamily = []
	Genus = []
	Species = []
	Subspecies = []
	Other = []

	query_phylogeny = []
	query_genus = ''
	query_name = ''

	# find the phylogeny of the query sequence
	# % identity = 100
	for hit in hit_accessions:

		if hit[1] == 100.0:

			query_phylogeny_string = hit[2] + ';'

			# parses the phylogeny string into a list (query_phylogeny)
			last_break = 0
			for index, c in enumerate(query_phylogeny_string):

				if c == ';':

					current_level = query_phylogeny_string[last_break:index]

					query_phylogeny.append(current_level)

					last_break = index + 1

			# add query name to the variable
			query_name = hit[3].splitlines()[0]

	# determine the species and subspecies strings for the query sequence
	#query_species, query_subspecies = get_deep_phylogeny(query_name)
	query_genus = query_phylogeny[-1]

	# find and compare phylogeny of current hit to phylogeny of query
	for hit in hit_accessions:

		hit_phylogeny = []
		hit_name = hit[3].splitlines()[0]
		
		hit[2] = hit[2] + ';'

		# parses the phylogeny string into a list (query_phylogeny)
		last_break = 0
		for index, c in enumerate(hit[2]):

			if c == ';':

				current_level = hit[2][last_break:index]

				hit_phylogeny.append(current_level)

				last_break = index + 1
				last_break = index + 1

		# --------------------------------------------------------------------
		# compares the hit phylogeny list to the query phylogeny list

		closest_relation = compare_hit_to_query_phylogeny(hit_phylogeny, query_phylogeny, hit_name, query_species, query_subspecies)

		if closest_relation == 'subspecies':

			Subspecies.append(hit)

		elif closest_relation == 'species':

			Species.append(hit)

		elif closest_relation == 'genus':

			Genus.append(hit)

		elif closest_relation == 'subfamily':

			Subfamily.append(hit)

		elif closest_relation == 'family':

			Family.append(hit)

		else:

			Other.append(hit)

	# ---------------------------------
	# output hit file data

	if input_type > 'accession':

		output_path = '{0}PHYLO_DATA/'.format(result_output_directory(os.getcwd()))

		# make dir if it doesn't already exist
		if not os.path.exists(output_path):

			os.makedirs(output_path)

	else:

		output_path = '{0}PHYLO_DATA/'.format(output_file_path)

		# make dir if it doesn't already exist
		if not os.path.exists(output_path):

			os.makedirs(output_path)

	output_phylo(Subspecies, 'Subspecies', output_path)
	output_phylo(Species, 'Species', output_path)
	output_phylo(Genus, 'Genus', output_path)
	output_phylo(Subfamily, 'Subfamily', output_path)
	output_phylo(Family, 'Family', output_path)
	output_phylo(Other, 'Other', output_path)

	summary_data = summarize_query(Other, Family, Subfamily, Genus, Species, Subspecies, query_name, query_database)
	summary_path = output_path + 'summary.txt'

	summary_save = open(summary_path, 'w')
	summary_save.write(summary_data)
	summary_save.close

		

# compares the phylogenetic information from the query and hit and returns the closest level of 
# phylogenetic relation between the two sequences
def compare_hit_to_query_phylogeny(hit_phylogeny, query_phylogeny, hit_name, query_species, query_subspecies):

	# allows us to get at least the species / subspecies values for hits / queries with
	# malformed phylo data (from GenBank)

	# remove leading space from genus field
	if hit_phylogeny[-1][0] == ' ':

		hit_phylogeny[-1] = hit_phylogeny[-1][1:]

	# remove leading space from genus field
	if query_phylogeny[-1][0] == ' ':

		query_phylogeny[-1] = query_phylogeny[-1][1:]

	# get genus from query name to see if it is in the hit name
	query_genus = query_phylogeny[-1]

	if len(query_phylogeny) == 1:
		
		query_genus = get_genus_from_name(query_phylogeny[-1])

	# if the genus, species, and subspecies names are the same in query
	# ex: gallus gallus gallus
	if query_subspecies.replace(' ', '').lower() == query_species.replace(' ', '').lower() == query_phylogeny[-1].replace(' ', '').lower() and query_subspecies.replace(' ', '').lower() in hit_name.lower():

		hit_count = hit_name.lower().count(query_subspecies.replace(' ', '').lower())

		if hit_count == 3:

			return 'subspecies'

		# if genus and species are the same, but not the subspecies
		# ex: Gallus gallus spadiceus
		elif hit_count == 2 and hit_phylogeny[-1].replace(' ', '').lower() == query_phylogeny[-1].replace(' ', '').lower():

			return 'species'

		elif hit_count == 1:

			return 'genus'

	# if the species and subspecies are the same in the query
	# ex: Anguilla bicolor bicolor
	elif query_subspecies.replace(' ', '').lower() == query_species.replace(' ', '').lower() and query_subspecies.replace(' ', '').lower() in hit_name.lower():

		hit_count = hit_name.lower().count(query_subspecies.replace(' ', '').lower())

		if hit_count == 2:

			return 'subspecies'

		elif hit_count == 1:

			return 'species'

		else:

			return 'genus'


	# if there is a subspecies in the query, and it is in the hit name, they must be subspecies related
	elif len(query_subspecies) > 0 and query_subspecies in hit_name:

		return 'subspecies'

	# if there is just the query species in the hit name, they must be delegated to OTHER
	elif query_species in hit_name:

		return 'species'

	# handle in case the query sequence has a phylogeny length of 1
	# GenBank return error
	# 
	# in this case, query phylogeny contains the Genus, species, and subspecies names,
	# so we can find genus relationship by looking at simply the name of the query
	elif hit_phylogeny[-1] in query_phylogeny[0]:

		return 'genus'

	# handle in case both the query and the hit 'genus' field are outputting the
	# Genus species (subspecies) format
	elif len(hit_phylogeny) == 1:

		if len(query_subspecies) > 0 and query_subspecies in hit_phylogeny[0]:

			return 'subspecies'

		elif query_species in hit_phylogeny[0]:

			return 'species'

		elif query_genus in hit_phylogeny[0]:

			return 'genus'

	# if the data is not malformed (length of phylogeny > 1)
	# check genus, which will always be the last level in the phylogeny
	elif len(query_phylogeny) > 1 and len(hit_phylogeny) > 1:

		# stripping leading spaces from query and hit phylogeny information
		query_phylogeny[-1] = query_phylogeny[-1].replace(' ', '')
		hit_phylogeny[-1] = hit_phylogeny[-1].replace(' ', '')

		# Handling outdated GenBank info (previously separate and now combined genuses)

		if query_phylogeny[-1] == 'Cervus' and hit_phylogeny[-1] == 'Przewalskium':

			return 'genus'

		elif query_phylogeny[-1] == 'Przewalskium' and hit_phylogeny[-1] == 'Cervus':

			return 'genus'

		elif query_phylogeny[-1] == hit_phylogeny[-1]:

			return 'genus'

		# Since phylo data from Genbank isn't labeled, we need to find out specific levels of
		# relation using a list of necessary Family, SubFamily, and Genus names for our organisms
		# of interest
		else:

			# check reversed list of query phylogeny
			# (highest level of relatedness first)
			for level in query_phylogeny[::-1]:

				# strip leading spaces
				current_level = level.replace(' ', '')

				# if the current query level is a subfamily, check if that value is in the hit
				# ... if so, they are subfamily related
				if current_level in Subfamilies:

					# --------------------------------------------------
					# Check for situation where family and subfamily have same name in query but
					# not in hit

					# Boolean (True if there are repeated subfamilies)
					repeated_subspecies = query_repeated_subspecies_not_hit(current_level, query_phylogeny, hit_phylogeny)

					if repeated_subspecies:

						return 'family'


					for hit_level in hit_phylogeny[::-1]:

						current_hit_level = hit_level.replace(' ', '')

						if current_hit_level == current_level:

							return 'subfamily'

			# check for family relatedness
			for level in query_phylogeny[::-1]:

				current_level = level.replace(' ', '')

				if current_level in Families:

					for hit_level in hit_phylogeny[::-1]:

						current_hit_level = hit_level.replace(' ', '')

						if current_hit_level == current_level:

							return 'family'

	else:

		return 'other'

# checks if the current level is repeated in the query phylogeny but not in the hit phylogeny
# MEANING: family relation only between query and hit phylogeny
def query_repeated_subspecies_not_hit(current_level, query_phylogeny, hit_phylogeny):

	query_count = 0
	hit_count = 0

	for level in query_phylogeny:

		if level.replace(' ', '') == current_level:

			query_count += 1

	for level in hit_phylogeny:

		if level.replace(' ', '') == current_level:

			hit_count += 1

	return query_count != hit_count

# parses the genus name from the query name
def get_genus_from_name(query_name):

	for index, c in enumerate(query_name):

		if c == ' ':

			if query_name[0] == ' ':

				return query_name[1:index]

			else:

				return query_name[:index]

# summarizes the query information
#	Phylogenetic level
#		% divergence from query | accession number
#			(descending)
def summarize_query(Other, Family, Subfamily, Genus, Species, Subspecies, query_name, query_database):

	summary = ''

	summary_header = 'Percent div to common anc.|    Accession #  |  Seq Length  |   Name    '
	summary_header += '\n ---------------------------------------------- \n'

	summary += summary_header
	summary += summarize_level(Other, 'Other', query_name, query_database)
	summary += summarize_level(Family, 'Family', query_name, query_database)
	summary += summarize_level(Subfamily, 'Subfamily', query_name, query_database)
	summary += summarize_level(Genus, 'Genus', query_name, query_database)
	summary += summarize_level(Species, 'Species', query_name, query_database)
	summary += summarize_level(Subspecies, 'Subspecies', query_name, query_database)

	# TODO consolidate excel outputs into single file
	#consolidate_result_files()



	return summary

def summarize_level(hits, level, query_name, query_database):

	level_summary = ''
	
	ordered_hits = [ ]

	# hits for excel output for data analysis
	excel_hits = [ ]

	# calculate the % difference as opposed to the % identity
	for hit in hits:

		percent_difference = (100 - hit[1]) / 2
		shortened_difference = ("%.2f" % percent_difference)
		accession = hit[0]

		# adds % difference, accession number, length of the sequence, and name of seq
		# to a list of ordered hits
		ordered_hits.append( [shortened_difference, accession, hit[3] , hit[3].splitlines()[0]] )

		# excel hits
		# Accession | Name | % difference | Seq length
		excel_hits.append( [accession, hit[3].splitlines()[0], shortened_difference, len(hit[3])] )

	# grossly hacky sort by field of list member
	# sorts members of the ordered hits by the % identity
	#	% identity is a string at this point of the float it was in percent difference
	#	casted to float and multiplied by -1 so that sorting results in the highest % first
	# 
	# result is a list ordered by the highest % difference
	ordered_hits.sort(key=lambda x: (-1 * float(x[0])))

	# DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	# DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG

	# output the hits to excel file
	excel_hits.sort(key=lambda x: (-1 * float(x[2])))

	book = xlwt.Workbook(encoding = 'utf-8')

	sheet1 = book.add_sheet('Sheet 1')

	# set the column labels
	sheet1.write(0, 0, 'Accession')
	sheet1.write(0, 1, 'Name')
	sheet1.write(0, 2, 'Percent_div_to_common_ancestor')
	sheet1.write(0, 3, 'Seq_length')
	sheet1.write(0, 4, 'Level')

	for index, hit in enumerate(excel_hits):

		sheet1.write(index + 1, 0, hit[0])
		sheet1.write(index + 1, 1, hit[1].replace(',', '')) # removes any commas that would disrupt csv parsing
		sheet1.write(index + 1, 2, hit[2])
		sheet1.write(index + 1, 3, hit[3])
		sheet1.write(index + 1, 4, level)

	save_path = ''

	if input_type == 'accession':

		save_path = '{0}Excel_Hits/'.format(result_output_directory(os.getcwd()))
		
		# make dir if it doesn't already exist
		if not os.path.exists(save_path):

			os.makedirs(save_path)

	# cross_blasting (file_path given)
	else:

		save_path = '{0}Excel_Hits/'.format(output_file_path)

		# make dir if it doesn't already exist
		if not os.path.exists(save_path):

			os.makedirs(save_path)

	# save the file

	file_save = '{0}/{1}_hits.csv'.format(save_path, level)
	book.save(file_save)

	# convert xls to csv
	with xlrd.open_workbook(file_save) as wb:

		sh = wb.sheet_by_index(0)
		
		with open(file_save, 'wb') as f:

			c = csv.writer(f)

			for r in range(sh.nrows):

				c.writerow(sh.row_values(r))

	# DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	# DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG

	# format to string
	level_summary += level + '\n'

	for hit in ordered_hits:

		# get the length of the FASTA seq (with newlines stripped)
		hit_seq = hit[2].splitlines()[1:]
		hit_seq = ''.join(hit_seq)
		hit_seq.replace('\n', '')
		seq_len = len(hit_seq)

		# value to add to summary file
		# % difference  | accession number | length of fasta seq
		line = '    {0}    |    {1}    |    {2}    |    {3}{4}'.format(hit[0], hit[1], seq_len, hit[3][28:], '\n')
		level_summary += line

	return level_summary

# outputs the final info for the identified hits
def output_phylo(hits, level, path):

	level_path = path + level + '/'

	# outputs the raw data for any processing that needs to be done
	for hit in hits:

		hit_path = '{0}{1}_{2}/'.format(level_path, hit[1], hit[0])

		# make dir if it doesn't already exist
		if not os.path.exists(hit_path):

			os.makedirs(hit_path)

		fasta_path = '{0}_fasta.txt'.format(hit_path)
		phylo_path = '{0}_phylo.txt'.format(hit_path)

		fasta_save = open(fasta_path, 'w')
		fasta_save.write('{0}'.format(hit[3]))
		fasta_save.close()

		# saves the phylo and % identity info
		phylo_save = open(phylo_path, 'w')
		phylo_save.write('{0}'.format(hit[1]))
		phylo_save.write('\n')
		phylo_save.write('{0}'.format(hit[2]))
		fasta_save.close()

# requires user interaction to determine the query sequence's species and subspecies strings
#		NOT available through GenBank query
def get_deep_phylogeny(query_name):

	querying_species = True
	querying_subspecies = True

	species = ''
	subspecies = ''

	while querying_species:

		print '\n Query Sequence: {0}'.format(query_name)
		current_species = raw_input('What is the species name of this organism? \n')
		correct = raw_input('Is this the correct species name? (y / n): {0} \n'.format(current_species))

		if correct.lower() == 'y' or correct.lower == 'yes':

			species = current_species

			querying_species = False

	while querying_subspecies:

		is_subspecies = raw_input('Is there a subspecies name of this organism? \n')

		if is_subspecies.lower() == 'y' or is_subspecies.lower() == 'yes':

			current_subspecies = raw_input('What is the subspecies name of this organism? \n')
			sub_correct = raw_input('Is this the correct subspecies name? (y / n): {0} \n'.format(current_subspecies))

			if sub_correct.lower() == 'y' or sub_correct.lower() == 'yes':

				subspecies = current_subspecies

				querying_subspecies = False

		else:

			querying_subspecies = False


	return species, subspecies




# queries GenBank for all of the hits in hit_accessions (with duplicates removed) to
# obtain the FASTA sequence and phylogenetic information for each hit
#		also returns the list of all accession numbers for file output
def query_genbank(query_name, query_database):

	querying = True

	# aggregator for all the accession numbers
	all_accessions = ''
	# aggregator for all of the fetched FASTA sequences
	full_hit_seqs = ''

	while querying:

		try:

			for index, hit in enumerate(hit_accessions):

				accession_num = hit[0]

				# saves the accession number in file for debugging or future reference
				all_accessions += (accession_num + '\n')

				# look up the accession number in GenBank database and download
				# the fasta sequence for that accession number
				handle = Entrez.efetch(db='nucleotide', id=accession_num, rettype='fasta')

				# ------------- PHYLOGENETIC INFORMTION -----------------
				# -------------------------------------------------------
				# get the organism phylogeny from genbank
				phylo_handle = Entrez.efetch(db='nucleotide', id=accession_num, rettype='gb', retmode='xml')

				if input_type == 'accession':

					# output file to folder "GENBANK_DATA", 1 level back from pwd
					file_path = '{0}GENBANK_DATA/'.format(result_output_directory(os.getcwd()))

					# make dir if it doesn't already exist
					if not os.path.exists(file_path):

						os.makedirs(file_path)

					file_name = '{0}query_{1}_hit_{2}.xml'.format(file_path, query_name, accession_num)

					save_file = open(file_name, 'w')
					save_file.write(phylo_handle.read())
					save_file.close()

					phylo_path = get_phylo(file_name, query_database)

					# get the full fasta sequence of the hit (not just the homologous region)
					current_fasta = handle.read()

					# add the fasta sequence and phylo path to hit_accessions
					hit_accessions[hit_accessions.index(hit)].append(phylo_path)
					hit_accessions[hit_accessions.index(hit)].append(current_fasta)

				# for cross_blasting
				else:

					file_path = '{0}GENBANK_DATA/'

					# make dir if it doesn't already exist
					if not os.path.exists(file_path):

						os.makedirs(file_path)

					file_name += 'query_{1}_hit_{2}.xml'.format(output_file_path, query_name, accession_num)

					save_file = open(file_name, 'w')
					save_file.write(phylo_handle.read())
					save_file.close()

					phylo_path = get_phylo(file_name, query_database)

					# get the full fasta sequence of the hit (not just the homologous region)
					current_fasta = handle.read()

					# add the fasta sequence and phylo path to hit_accessions
					hit_accessions[hit_accessions.index(hit)].append(phylo_path)
					hit_accessions[hit_accessions.index(hit)].append(current_fasta)

				# -------------------------------------------------------
				# -------------------------------------------------------

				

				full_hit_seqs += current_fasta

			return all_accessions, full_hit_seqs

		except:

			query_genbank(query_name, query_database)

		querying = False

# consolidates duplicate accession numbers in the BLAST query hits to get accurate identity %s
def remove_duplicate_hits(accessions):

	previous_numbers = {}

	hit_accession_replacements = []

	for hit in accessions:

		if not hit[0] in previous_numbers:

			previous_numbers[hit[0]] = hit[1]

	x = 0

	for accession, identity in previous_numbers.iteritems():

		hit_accession_replacements.append([accession, identity])

		x += 1

	# trim any remaining / duplicate values from hit_accessions
	hit_accessions[:] = hit_accession_replacements
	

# returns the phylogenetic information of the given hit
def get_phylo(xml_file, query_database):

	# grabs the phylogeny of the hit sequence
	tree = ET.parse(xml_file)
	root = tree.getroot()

	# location is different if seq is a reference genome
	if query_database == 'refseq':

		phylo = root[0][16].text

	else:

		phylo = root[0][14].text

	return phylo

# returns just the GenBank accession number from a FASTA header
def get_accession(title):

	first_index = 0
	final_index = 0

	bar_count = 0

	for index, c in enumerate(title):		

		if c == '|':

			bar_count += 1

		if bar_count == 2:

			first_index = index + 2

		if bar_count == 3:

			final_index = index + 1

	# strip the accession from the title
	return title[first_index:final_index]

# outputs a single file of the hits FASTA sequences
def output_full_hits(seq_name, all_hits, all_accessions):

	# prevent IO string length error from improperly named GenBank entries
	if len(seq_name) > 75:

		seq_name = seq_name[:75]

	if input_type == 'accession':
		
		file_path = result_output_directory(os.getcwd())

		# make dir if it doesn't already exist
		if not os.path.exists(file_path):

			os.makedirs(file_path)

		result_file = '{0}{1}_FASTA.txt'.format(file_path, 'full_hits')

		print '\nWriting output files ...'

		save_file = open(result_file, 'w')
		save_file.write(all_hits)
		save_file.close()

	else:

		result_file = '{0}{1}_FASTA.txt'.format(output_file_path, 'full_hits')

		print '\nWriting output files ...'

		save_file = open(result_file, 'w')
		save_file.write(all_hits)
		save_file.close()

# strips a leading unicode defining character 'u' from the hit data
def strip_incorrect_chars(hit):

	new_hit = []

	for i in hit:

		i = i[1:]

		new_hit.append(i)

	return new_hit

# returns the pwd, minus three levels of depth
def result_output_directory(current_directory):

	rev_dir = current_directory[::-1]

	rev_result = ''

	result = ''

	count = 0

	if input_type == 'accession':

		for index, c in enumerate(rev_dir):

			if count == 2:

				rev_index = len(current_directory) - (index)
				
				result = current_directory[:rev_index]

				result += '/Results/{0}/SingleBLAST_{1}_{2}_({3}h_{4}m_{5}s)/{6}_{7}_{8}/'.format(today, query_species, query_subspecies, hour, minute, second, query_species, query_subspecies, query_database)

				return result

			elif c == '/':

				count += 1

# TODO
# consolidates all of the Excel output files into a single summary output csv file
def consolidate_result_files():

	root = Tkinter.Tk()
	root.withdraw()
	file_paths = tkFileDialog.askopenfilenames()

	for index, path in enumerate(file_paths):

		print '' 

# runs the business
def main():

	# sets values to global for naming dirs for result output
	global input_type
	global query_species
	global query_subspecies
	global query_database
	global query_accession
	global output_file_path

	# defines the type of run that will be completed
	# i.e. FASTA, accession #, or w.e.
	input_type = sys.argv[1]
	query_database = sys.argv[2]
	query_species = sys.argv[3]
	query_subspecies = sys.argv[4]
	query_accession = sys.argv[5]
	output_file_path = sys.argv[6]

	# creates the result directory if it doesn't exist
	if not os.path.exists(output_file_path):

		os.makedirs(output_file_path)

	if input_type == 'F' or input_type == 'f':

		print '\n----- FASTA Input -----\n '

		select_files()
		handle_requests(query_database, 'fasta', None)
		parse_queries(query_database, query_species, query_subspecies, input_type)

		print '\n----- Complete! -----\n '

	elif input_type == 'A' or input_type == 'a':

		print '\n---- Accession Input -----\n '

		handle_requests(query_database, 'accession', query_accession)
		parse_queries(query_database, query_species, query_subspecies, input_type)

		print '\n----- Complete! -----\n '

	elif input_type == 'cross':

		print '\n---- Cross BLAST Accession Input -----\n '

		handle_requests(query_database, 'cross', query_accession)
		parse_queries(query_database, query_species, query_subspecies, input_type)

		print '\n----- Complete! -----\n '

	elif input_type == 'genbank':

		print '\n---- GenBank Query with BLAST XML Input -----\n '

		root = Tkinter.Tk()
		root.withdraw()
		file_path = tkFileDialog.askopenfilename()
		xml_results.append(file_path)
		parse_queries(query_database, query_species, query_subspecies, input_type)

		print '\n----- Complete! -----\n '

	else:

		raise Exception("----- ERROR: Valid input type needed -----")

if __name__ == "__main__":

    main()

'''

TODO

	[ ] ensure that all is working after updates to cross_blast -> blast_accession pipeline (esp. with file output)
		[ ] need to make directories when they don't exist (inside the result dir)
	[ ] Set global variables / shell arguments to default values that get overridden if user inputs values
	[ ] Add species and subspecies data to csv summary files
	[ ] Make so that sequences that only have species and no subspecies match on genus level, not species level
		- may need to make a part of manual data parsing with result .csv files
	[ ] Figure out how to re-enter species / subspecies names if mistyped that doesn't
		result in mis-formatted file-output and query_species / query_subspecies data
		entries.

'''


