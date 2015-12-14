import Tkinter, tkFileDialog
import csv
import os
import datetime

# date and time information for writing to file
today = datetime.datetime.today().strftime('%y_%m_%d')
hour = datetime.datetime.today().time().hour
minute = datetime.datetime.today().time().minute
second = datetime.datetime.today().time().second

class Condenser:

	def __init__(self, files):

		self.files = files

	# selects a number of files that will be condensed into a single file
	def collect_files(self):

		root = Tkinter.Tk()
		root.withdraw()
		file_paths = tkFileDialog.askopenfilenames()

		for file in file_paths:

			self.files.append(file)

	# condenses the contents of the selected files to a single file
	def condense_results(self):

		dir_path = one_directory_back(os.getcwd())
		dir_path += 'Summaries/'

		# make dir if it doesn't already exist
		if not os.path.exists(dir_path):

			os.makedirs(dir_path)

		save_file = '{0}Summary_{1}_({2}h_{3}m_{4}s).csv'.format(dir_path, today, hour, minute, second)

		with open(save_file, 'wb') as f:

			c = csv.writer(f, delimiter = ',', quotechar = '|')

			for index, csv_file in enumerate(self.files):

				if index == 0:

					with open(csv_file) as csvfile:

						reader = csv.reader(csvfile, delimiter = ',', quotechar = '|')

						for row in reader:

							c.writerow(row)

				else:

					with open(csv_file) as csvfile:

						reader = csv.reader(csvfile, delimiter = ',', quotechar = '|')

						for index, row in enumerate(reader):

							if index > 0:
								
								c.writerow(row)


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

def main():

	os.system('clear')

	condenser = Condenser([])
	condenser.collect_files()
	condenser.condense_results()

if __name__ == '__main__':

	main()















