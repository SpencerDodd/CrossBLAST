import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

test_file = '/Users/spencerdodd/Desktop/Summary_15_11_28_(14h_18m_10s).csv'

class HistogramParser:

	def __init__(self, other, family, subfamily, genus, species, subspecies):

		self.other = other
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

					phylo_level = row[5]
					percent_div_to_common_anc = float(row[3])

					if phylo_level == 'Other':

						self.other.append(percent_div_to_common_anc)

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

		level_names = ['Other', 'Family', 'Subfamily', 'Genus', 'Species', 'Subspecies']
		levels = [self.other, self.family, self.subfamily, self.genus, self.species, self.subspecies]
		colors = ['dodgerblue', 'cornflowerblue', 'deepskyblue', 'turquoise', 'cyan', 'lightgreen']
		bin_max = 10

		# finds the max value in the data set to set histogram x-range size
		for level in levels:

			if len(level) > 0:

				level_max = max(level)

				if level_max > bin_max:

					bin_max = level_max

		bins = np.linspace(0, bin_max, 100)

		for index, level in enumerate(levels):

			if len(level) > 0:
				
				plt.hist(level, bins, alpha = 0.5, label = level_names[index])

		plt.ylabel('Frequency')
		plt.xlabel('Percent dist to common ancestor')
		plt.title('Overview')
		plt.legend(loc = 'upper right')
		plt.savefig('/Users/spencerdodd/Desktop/Overview.png')

def main():

	parser = HistogramParser([], [], [], [], [], [])
	parser.hist_parse(test_file)
	parser.make_hists()

main()





