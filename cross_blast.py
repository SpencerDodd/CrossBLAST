import sys
import os
import Tkinter, tkFileDialog
import csv
from subprocess import call
from Bio import Entrez

# blasts all species-related sequences and accumulates distances for calulation
# 	of minimums and maximums

class CrossBlast:

	def __init__(self, request_path, accessions, query_names):

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

		for index, accession in enumerate(self.accessions):

			print 'Sequence {0} of {1}'.format(index + 1, len(self.accessions))

			handle = Entrez.efetch(db='nucleotide', id=self.accessions[index], rettype='fasta')

			record = handle.read()

			name = record.splitlines()[0]

			query_species, query_subspecies = get_deep_phylogeny(name)

			self.query_names.append( [query_species, query_subspecies] )

	# blasts all of the collected accessions
	def blast_accessions(self):

		for index, accession in enumerate(self.accessions):

			accession = self.accessions[index]
			species = self.query_names[index][0]
			subspecies = self.query_names[index][1]

			call(['python', 'BLAST_Accession.py', 'cross', species, subspecies, accession])





	# BLASTs all of the accessions 

# returns the pwd, minus one level of depth
def one_directory_back(current_directory):

	rev_dir = current_directory[::-1]

	rev_result = ''

	result = ''

	for c in rev_dir:

		if c == '/':

			rev_result += rev_dir[rev_dir.index(c):]
			
			result = rev_result[::-1]

			return result

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
	
def main():

	# clears the current terminal shell
	os.system("clear")

	request = CrossBlast(None, [], [])
	request.select_sequences()
	request.get_accessions()
	request.get_phylo_info()
	request.blast_accessions()

if __name__ == '__main__':

	main()


# TODO
#	[ ] Change BLAST_Accession.py to input species and subspecies info at the beginning
#		of the blasting as a main[x] input








