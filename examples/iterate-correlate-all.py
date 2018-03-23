""" iterate-correlate-all.py: Example usage of the rpkm_correlate library - iterates through all PWY/ORF files in two directories and correlates them, putting
							  the output in a third directory. 

							  Requirements:
							  - rpkm_correlate.py in the current working directory 
							  - data/ folder in the current working directory
							  - pwy/ folder within the data folder containing pathway info files named as maxbin_###.pwy.txt
							  - orf/ folder within the data folder containing ORF RPKM data files named as maxbin_###.orf_rpkm.txt
							  - out/ folder within the data folder
							  - referenced files are properly formatted (see README.md)

							  """

__author__ = "Eli Kaplan"
__email__ = "eli.kaplan@alumni.ubc.ca"
__license__ = "CC0 - No rights reserved. See LICENSE"

from rpkm_correlate import correlateRPKM

# os.listdir and os.path.join are necessary to generate a list of files to operate on
from os import listdir
from os.path import join

# Generate the list of pathway info files to work with
pwy_files = [join('data/pwy/', f) for f in listdir('data/pwy/')]


# Iterate through and correlate all of the relevant files
for file in pwy_files:
	# Generate the name of the ORF data file that corresponds to the input PWY file
	orf_name = file.replace('.pwy.txt', '.orf_rpkm.txt').replace('data/pwy/','data/orf/')

	# Generate a fitting name for the output file
	out_name = file.replace('data/pwy/', 'data/out/').replace('.pwy.txt','_out.txt')

	print(file + ' + ' + orf_name + ' -> ' + out_name)

	# Correlate the data and output it. 
	correlateRPKM(file, orf_name, out_name)
