import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

class HistogramParser:

	def __init__(self, other, superfamily, family, subfamily, genus, species, subspecies):

		self.other = other
		self.superfamily = superfamily
		self.family = family
		self.subfamily = subfamily
		self.genus = genus
		self.species = species
		self.subspecies = subspecies

	# parses the input file results by phylogenetic group
	def hist_parse(self, input_file):

		with open(input_file) as csv_file:

			reader = csv.reader(csv_file, delimiter = ',', quotechar = '|')

			for index, row in enumerate(reader):

				if index != 0:

					phylo_level = row[4]
					percent_div_to_common_anc = float(row[2])

					if phylo_level == 'Other':

						self.other.append(percent_div_to_common_anc)

					elif phylo_level == 'Superfamily':

						self.superfamily.append(percent_div_to_common_anc)

					elif phylo_level == 'Family':

						self.family.append(percent_div_to_common_anc)

					elif phylo_level == 'Subfamily':

						self.subfamily.append(percent_div_to_common_anc)

					elif phylo_level == 'Genus':

						self.genus.append(percent_div_to_common_anc)

					elif phylo_level == 'Species':

						self.species.append(percent_div_to_common_anc)

					elif phylo_level == 'Subspecies':

						self.subspecies.append(percent_div_to_common_anc)

					else:

						raise Exception('ERROR | Phylogenetic level of hit not properly defined: {0}'.format(phylo_level))

	# creates and saves the histograms for all phylogenetic levels of the results summary 
	def make_hists(self):

		level_names = ['Other', 'Superfamily', 'Family', 'Subfamily', 'Genus', 'Species', 'Subspecies']
		levels = [self.other, self.superfamily, self.family, self.subfamily, self.genus, self.species, self.subspecies]
		colors = ['dodgerblue', 'blue','cornflowerblue', 'deepskyblue', 'turquoise', 'cyan', 'lightgreen']
		bin_max = 10

		# finds the max value in the data set to set histogram x-range size
		for level in levels:

			if len(level) > 0:

				level_max = max(level)

				if level_max > bin_max:

					bin_max = level_max

		bins = np.linspace(0, bin_max, 200)

		for index, level in enumerate(levels):

			if len(level) > 0:
				
				plt.hist(level, bins, alpha = 0.5, label = '{0} (Range: {1} to {2})'.format(level_names[index], min(level), max(level)))

		plt.ylabel('Frequency')
		plt.xlabel('Percent dist to common ancestor')
		plt.title('Overview')
		plt.legend(loc = 'upper right')
		plt.savefig('{0}Overview.png'.format(results_dir))

def main():

	global results_dir
	results_dir = sys.argv[1]

	parser = HistogramParser([], [], [], [], [], [], [])
	parser.hist_parse('{0}Summary.csv'.format(results_dir))
	parser.make_hists()

main()

'''
TODO

'''








