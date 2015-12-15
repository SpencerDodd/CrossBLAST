import sys
import os
import Tkinter, tkFileDialog
import csv
from subprocess import call
from Bio import Entrez
import datetime
from twilio.rest import TwilioRestClient

start_time = datetime.datetime.now()
today = datetime.datetime.today().strftime('%y_%m_%d')
hour = datetime.datetime.today().time().hour
minute = datetime.datetime.today().time().minute
second = datetime.datetime.today().time().second

# blasts all species-related sequences and accumulates distances for calulation
# 	of minimums and maximums

class CrossBlast:

	def __init__(self, query_database, request_path, accessions, query_names):

		self.query_database = query_database
		self.request_path = request_path
		self.accessions = accessions
		self.query_names = query_names

	# selects the output file with the sequences to blast
	def gather_sequences(self):

		print '\nGathering initial BLAST results ...\n'

		self.request_path = initial_file_path + 'Excel_hits/Species_hits.csv'

	# the initial BLAST run to generate the Species.csv that is used for cross_blasting
	def initial_blast(self):
		
		handle = Entrez.efetch(db='nucleotide', id=initial_accession, rettype='fasta')
		record = handle.read()
		name = record.splitlines()[0]

		global initial_species, initial_subspecies
		initial_species, initial_subspecies = get_deep_phylogeny(name)

		global initial_file_path
		initial_file_path = output_initial_directory()

		print '\n---- Querying Initial BLAST Sequence: {0} {1} ({2}) ----'.format(initial_species, initial_subspecies, initial_accession)
		
		# update with progress
		call(['python', 'send_update.py', 'CrossBLAST (initial sequence) {0}: {1} {2}'.format(initial_accession, initial_species, initial_subspecies)])
		# BLASTs
		call(['python', 'blast_accession.py', 'cross', query_database, initial_species, initial_subspecies, initial_accession, initial_file_path])


	# reads the accession numbers from the results file and adds to the accessions
	# parameter of the CrossBlast object
	def get_accessions(self):

		with open(self.request_path, 'rb') as csvfile:

			data = csv.reader(csvfile, delimiter = ',')

			for index, row in enumerate(data):

				if index > 0:

					print row[0]
					self.accessions.append(row[0])

	# gets the species and subspecies names of each accession
	def get_phylo_info(self):

		query_species_name = ''

		for accession in self.accessions:
			
			print 'Sequence {0} of {1}'.format(index + 1, len(self.accessions))

			handle = Entrez.efetch(db='nucleotide', id=self.accessions[index], rettype='fasta')
			record = handle.read()
			name = record.splitlines()[0]
			query_subspecies = get_deep_phylogeny(name)

			self.query_names.append( [query_species_name, query_subspecies] )

	# blasts all of the collected accessions
	def blast_accessions(self):

		for index, accession in enumerate(self.accessions):

			accession = self.accessions[index]
			species = self.query_names[index][0]
			subspecies = self.query_names[index][1]
			file_path = output_directory(species, subspecies, query_database)

			# clears the current terminal shell
			os.system("clear")

			print 'Querying sequence {0} / {1}'.format(index + 1, len(self.accessions))

			# update with progress
			call(['python', 'send_update.py', 'CrossBLAST-ing {0}: {1} {2} (Sequence {3} / {4}'.format(accession, species, subspecies, index + 1, len(self.accessions))])
			# BLAST
			call(['python', 'blast_accession.py', 'cross', self.query_database, species, subspecies, accession, file_path])

# returns the output directory of the initial query
def output_initial_directory():

	current_directory = os.getcwd()

	rev_dir = current_directory[::-1]

	rev_result = ''

	result = ''

	count = 0

	for index, c in enumerate(rev_dir):

		if count == 2:

			rev_index = len(current_directory) - (index)
			
			result = current_directory[:rev_index]

			# variables for filenaming
			cross_blast_query_name = str(initial_species) + '_' + str(initial_subspecies)
			
			result += '/Results/{0}/CrossBLAST_{1}_({2}h_{3}m_{4}s)/{5}/'.format(today, cross_blast_query_name, hour, minute, second, 'Initial')

			return result

		elif c == '/':

			count += 1

# returns the pwd, minus three levels of depth
def output_directory(species, subspecies, query_database):

	current_directory = os.getcwd()

	rev_dir = current_directory[::-1]

	rev_result = ''

	# so that we can use this dir to condense the results and produce histograms from
	global result
	result = ''

	count = 0

	for index, c in enumerate(rev_dir):

		if count == 2:

			rev_index = len(current_directory) - (index)
			
			result = current_directory[:rev_index]

			# variables for filenaming
			cross_blast_query_name = str(species) + '_' + str(subspecies)
			
			result += '/Results/{0}/CrossBLAST_{1}_({2}h_{3}m_{4}s)/{5}_{6}_{7}/'.format(today, cross_blast_query_name, hour, minute, second, species, subspecies, query_database)

			return result

		elif c == '/':

			count += 1

# auto-parses phylo information from the given string
# returns second two words separated by spaces:
# Example:
#		given: >gi|62184368|ref|NC_006915.1| Mus musculus molossinus mitochondrion, complete genome
# 		returns: musculus, molossinus
def get_deep_phylogeny(query_name):

	phylogeny = []

	column_count = 0
	space_count = 0

	for index, c in enumerate(query_name):

		if column_count == 4:

			query_remainder = query_name[index + 1:]

			break_point = 0

			for index, c in enumerate(query_remainder):

				if space_count == 3:

					# returns the 2nd and 3rd words in the name
					# species and subspecies (or garbage for subspecies)
					return phylogeny[1], phylogeny[2]
					break

				elif c == ' ':

					phylo_section = query_remainder[break_point:index]
					phylogeny.append(phylo_section)
					break_point = index + 1
					space_count += 1

		elif c == '|':

			column_count += 1

def main():

	# clears the current terminal shell
	os.system("clear")

	global query_database
	global initial_accession
	query_database = 'refseq_genomic'
	initial_accession = sys.argv[1]

	request = CrossBlast(query_database, None, [], [])

	# initial BLAST run to get the species.csv
	request.initial_blast()
	# now run the rest of the CrossBLAST
	request.gather_sequences()
	request.get_accessions()
	request.get_phylo_info()
	request.blast_accessions()

	# condenses the result files automatically into a single csv
	condense_results(output_directory(initial_species, initial_subspecies, query_database))

	# creates histogram image for the result
	hist_results(output_directory(initial_species, initial_subspecies, query_database))

	end_time = datetime.datetime.now()

	update_message = 'RUN COMPLETE! {0}: {1} {2}'.format(initial_accession, initial_species, initial_subspecies)
	update_message += '\n-----------'
	update_message += '\nRuntime: {0}'.format(end_time - start_time)
	call(['python', 'send_update.py', update_message])

if __name__ == '__main__':

	main()

'''

TODO:

	[ ] add a log file so that if you leave it running on a box and SSH in you can
		see the progress on the run / what remains

		- could also update the sequence's BLAST folder name with a (COMPLETED) tag



'''








