import sys
import os
import Tkinter, tkFileDialog
import csv
from subprocess import call
from Bio import Entrez
import datetime

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
	def select_sequences(self):

		print '\nSelect the .csv output file with the accessions you would like to blast\n'

		root = Tkinter.Tk()
		root.withdraw()
		file_path = tkFileDialog.askopenfilename()

		self.request_path = file_path

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

		for index, accession in enumerate(self.accessions):

			if index == 0:

				print 'Sequence {0} of {1}'.format(index + 1, len(self.accessions))

				handle = Entrez.efetch(db='nucleotide', id=self.accessions[index], rettype='fasta')

				record = handle.read()

				name = record.splitlines()[0]

				query_species_name, query_subspecies = get_deep_phylogeny(name, index)

				self.query_names.append( [query_species_name, query_subspecies] )

			else:

				print 'Sequence {0} of {1}'.format(index + 1, len(self.accessions))

				handle = Entrez.efetch(db='nucleotide', id=self.accessions[index], rettype='fasta')

				record = handle.read()

				name = record.splitlines()[0]

				query_subspecies = get_subspecies(name)

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

			call(['python', 'blast_accession.py', 'cross', self.query_database, species, subspecies, accession, file_path])

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

# requires user interaction to determine the query sequence's species and subspecies strings
#		NOT available through GenBank query
def get_deep_phylogeny(query_name, index):

	global origin_species
	global origin_subspecies

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
			origin_species = current_species

			querying_species = False

	while querying_subspecies:

		is_subspecies = raw_input('Is there a subspecies name of this organism? \n')

		if is_subspecies.lower() == 'y' or is_subspecies.lower() == 'yes':

			current_subspecies = raw_input('What is the subspecies name of this organism? \n')
			sub_correct = raw_input('Is this the correct subspecies name? (y / n): {0} \n'.format(current_subspecies))

			if sub_correct.lower() == 'y' or sub_correct.lower() == 'yes':

				subspecies = current_subspecies
				origin_subspecies = current_subspecies

				querying_subspecies = False

		else:

			querying_subspecies = False


	return species, subspecies

# just asks for the subspecies name if it isn't the first sequence
def get_subspecies(query_name):

	querying_subspecies = True

	subspecies = ''

	while querying_subspecies:

		print '\n Query Sequence: {0}'.format(query_name)
		current_subspecies = raw_input('What is the subspecies name of this organism? \n')
		sub_correct = raw_input('Is this the correct subspecies name? (y / n): {0} \n'.format(current_subspecies))

		if sub_correct.lower() == 'y' or sub_correct.lower() == 'yes':

			subspecies = current_subspecies

			querying_subspecies = False

		else:

			querying_subspecies = False


	return subspecies

	
def main():

	# clears the current terminal shell
	os.system("clear")

	global query_database
	query_database = sys.argv[1]

	request = CrossBlast(query_database, None, [], [])
	request.select_sequences()
	request.get_accessions()
	request.get_phylo_info()
	request.blast_accessions()

	# condenses the result files automatically into a single csv
	condense_results(output_directory(species, subspecies, query_database))

	# creates histogram image for the result
	hist_results(output_directory(species, subspecies, query_database))

if __name__ == '__main__':

	main()

'''

TODO

	[ ] species will be the same...don't ask for it, just subspecies
	[ ] add 'y/n' options to presence of subspecies
	[ ] make so that full cross_blast will run off of a single accession number
	[ ] tie condense_results.py and hist_results.py into this script



'''








