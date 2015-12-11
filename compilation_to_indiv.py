# Script to output individual seq files from a compilation of FASTA or MAFFT seqs
# TODO code needs to be mutated
'''
import Tkinter, tkFileDialog

def remove_alignment():

	root = Tkinter.Tk()
	root.withdraw()
	file_paths = tkFileDialog.askopenfilenames()

	for path in root.tk.splitlist(file_paths):

		with open (path, 'r') as myfile:

			final_data = ''

			data = myfile.read().splitlines(True)

			for l in data:

				l = l.replace('-', '')

				final_data = final_data + l

		file_name = path + 'FASTA.txt'

		final_data = reformat_text(final_data)

		out_file = open(file_name, 'w')
		out_file.write(final_data)
		out_file.close()


def reformat_text(final_data):

	seqs = []

	return_data = ''

	current_seq = 0

	final_data = final_data.splitlines(True)

	for l in final_data[1:]:

		index = final_data.index(l)

		if len(l) > 0:

			if l[0] == '>':

				completed_seq = final_data[current_seq:index]
				seqs.append(completed_seq)
				
				current_seq = index

	for s in seqs:

		seq_string = ''

		for l in s:

			if l[0] != '>':
				
				l = l.strip('\r\n')

			if len(l) > 0:

				seq_string += l

		seq_string += '\r\n'

		return_data += seq_string

	return return_data




def main():

	remove_alignment()

main() '''