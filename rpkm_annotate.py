""" rpkm_annotate.py : Library / script to combine pathway information, RPKM data, and annotation information into output data"""

__author__ = "Eli Kaplan"
__email__ = "eli.kaplan@alumni.ubc.ca"
__license__ = "CC0 - No rights reserved. See LICENSE"

import csv
import sys

from os import listdir
from os.path import join

from rpkm_correlate import loadPathwayInfoFromFile, loadORFDataFromFile, correlatePathwayInfoWithData


def loadAnnotationsFromFile(filename, sample_name, csv_separator):
	""" loadAnnotationsFromFile(): 	Loads ORF annotations from a specified file and returns the relevant fields for each ORF

									Parameters:
									- filename: name of the file to load annotations from
									- sample_name: name of the sample these annotations are for (to properly format the IDs)
									- csv_separator: column separator used in given file

									Returns: a tuple - (name of sample, dict of tuples of annotation data indexed by ORF ID) """

	anno_data=(sample_name, {})

	try:
		# Load the specified annotations file
		with open(filename, 'r') as anno_file:
			anno_reader = csv.reader(anno_file, delimiter=csv_separator)

			# Process each row in the annotations file and store the relevant data
			for row in anno_reader:
				if row[0] == '#query': # Skip the header row
					continue

				anno_query = 'O_' + row[0].replace(sample_name,'')[1:] # Generate the proper corresponding ORF ID

				# Only take the first annotation result from each ORF
				if anno_query in anno_data[1]:
					continue

				anno_hit = row[9].split('[')[0] # Grab only the name of the gene as the 'hit'
				anno_q_length = row[2] 	# 'q_length' column
				anno_bitscore = row[3] 	# 'bitscore' column
				anno_bsr = row[4] 		# 'bsr' column
				anno_expect = row[5]	# 'expect' column
				anno_identity = row[7] 	# 'identity' column
				anno_ec = row[8] 		# 'ec' column

				# Store the loaded data in the dict in anno_data, which is keyed by the ORF ID
				anno_data[1][anno_query] = (anno_hit, anno_q_length, anno_bitscore, anno_bsr, anno_expect, anno_identity, anno_ec)

	# If an error was encountered loading and processing the annotation file, exit and print the error
	except Exception as e:
		print("Error reading annotation file - exiting.")
		print("Exception: " + str(e))
		quit()


	return anno_data


def batchCorrelateAnnotate(file_dir, output_filename = 'pwy_anno.tsv', csv_separator='\t', pwy_file_suffix='.pwy.txt', data_file_suffix='.orf_rpkm.txt', anno_file_suffix='.metacyc-2016-10-31.lastout.parsed.txt', selected_pathways=[]):
	""" batchCorrelateAnnotate() : 	Load pathway information, RPKM data, and annotation data from separate files in a given directory, and produce an output file showing the relevant annotations and RPKM value
									for each ORF in each pathway in each sample in the loaded files.

									Parameters:
									- file_dir: the path to the directory in which all of the relevant files are stored
									- output_filename: the name of the file to store the output of this function in
									- csv_separator: column separator used in input files and output file
									- pwy_file_suffix: suffix for files containing pathway information
									- data_file_suffix: suffix for files containing RPKM data
									- anno_file_suffix: suffix for files containing annotation data
									- selected_pathways: list of strings of pathway short names, specifiying which pathways should be included in the output file

									Returns: nothing - outputs results to specified file."""

	# Create a list of all files in the target directory
	all_files = [join(file_dir, f) for f in listdir(file_dir)]

	# Split the list of all the files into three lists: pathway files, data files, and annotation files
	pwy_files = []
	data_files = []
	anno_files = []

	for file in all_files:
		if pwy_file_suffix in file:
			pwy_files.append(file)
		elif data_file_suffix in file:
			data_files.append(file)
		elif anno_file_suffix in file:
			anno_files.append(file)
		#else:
			#print("Unknown file in batch directory: " + file)


	# Match each pathway information file with its corresponding RPKM data file and annotation file
	file_pairs = []

	for pwy_file in pwy_files:
		# Generate the name of the corresponding data file based on the name of the pathway info file
		corresponding_data_file = pwy_file.replace(pwy_file_suffix, data_file_suffix)

		# Generate the name of the corresponding annotation file, using the same technique
		corresponding_anno_file = pwy_file.replace(pwy_file_suffix, anno_file_suffix)
	
		# If the necessary files exist, pair them with the pathway file. Otherwise, don't process anything for this set of pathway data.
		if corresponding_data_file in data_files:
			if corresponding_anno_file in anno_files:
				file_pairs.append((pwy_file, corresponding_data_file, corresponding_anno_file))

			else:
				print("Missing annotation file: " + corresponding_anno_file)

		else:
			print("Missing data file: " + corresponding_data_file)


	# Produce the output data for this function
	output_data = []

	# Keep track of the total number of pathway/sample pairs, RPKM data points, and annotations loaded
	n_total_pathways = 0
	n_total_datapoints = 0
	n_total_annotations = 0

	for (pathway_file, data_file, anno_file) in file_pairs:
		# Load pathway information from the given file
		pathway_info = loadPathwayInfoFromFile(pathway_file, csv_separator)
		cur_sample = pathway_info[0]

		# Load RPKM data from the file corresponding to the pathway info file
		rpkm_data = loadORFDataFromFile(data_file, cur_sample, csv_separator)

		# Load annotation data from the corresponding file
		anno_data = loadAnnotationsFromFile(anno_file, cur_sample, csv_separator)

		# Create a dict from the RPKM data, indexed by ORF ID, for easier searching
		data_dict = dict(rpkm_data[1])

		# Keep track of the number of missing annotations and RPKM data points, for troubleshooting.
		n_missing_annotations = 0
		n_missing_rpkm = 0

		# Iterate through each pathway
		for pwy, pwy_cname, pwy_orfs in pathway_info[1]:

			# If a list of selected pathways is selected and this pathway is not in that list, skip it.
			if len(selected_pathways) > 0 and pwy not in selected_pathways:
				continue

			n_total_pathways += 1

			# Iterate through each ORF in this pathway
			for orf in pwy_orfs:

				# If the ORF is present in the RPKM data loaded for this sample, continue processing this ORF
				if orf in data_dict:
					n_total_datapoints += 1

					# If the annotation data is present for this ORF, add a new row for this ORF/pathway/sample pair to the output data
					if orf in anno_data[1]:
						n_total_annotations += 1

						orf_anno = anno_data[1][orf] # Acquire the annotation data for this specific ORF

						# Add a new row to the output data
						output_data.append((cur_sample, # SAMPLE
							pwy,			# PWY_NAME
							orf,			# ORF
							orf_anno[0], 	# HIT
							data_dict[orf], # RPKM
							orf_anno[1], 	# Q_LENGTH
							orf_anno[2],	# BITSCORE
							orf_anno[3],	# BSR
							orf_anno[4],	# EXPECT
							orf_anno[5],	# IDENTITY
							orf_anno[6]))	# EC


					else: # If the ORF is missing annotation data, output it but indicate that there are no annotations (with zeroes in values etc.)
						n_missing_annotations += 1

						output_data.append((cur_sample,
							pwy,
							orf,
							orf,
							data_dict[orf],
							0,
							0,
							0,
							0,
							0,
							0))


				else: # Keep track of missing RPKM data points
					n_missing_rpkm += 1

		print('Loaded sample: ' + cur_sample + ' - ORFS with no annotations: ' + str(n_missing_annotations) + ' - missing rpkm data points: ' + str(n_missing_rpkm))
				

	print('Processed ' + str(n_total_pathways) + ' sample/pathway pairs, ' + str(n_total_datapoints) + ' RPKM data points, ' + str(n_total_annotations) + ' total annotations.')


	
	# Generate a header for the output tabulated file 
	output_file_header = ['SAMPLE', 'PWY_NAME', 'ORF', 'HIT', 'RPKM', 'Q_LENGTH', 'BITSCORE', 'BSR', 'EXPECT', 'IDENTITY', 'EC']


	# Output the resulting data to the chosen file
	try:
		with open(output_filename, 'w') as output_file:
			output_writer = csv.writer(output_file, delimiter=csv_separator)

			# Write the header to the output file as the first line.
			output_writer.writerow(output_file_header)

			# Write out each row
			for row in output_data:
				output_writer.writerow(row)

	
	except Exception as e: # If an error occurred while writing out the results, exit.
		print("ERORR: File output failed - exiting.")
		print("Exception: " + str(e))
		quit()



# Allow the script to be run stand-alone (and prevent the following code from running when imported)
if __name__ == "__main__":
	def printUsage():
		"""Prints usage information for this script."""
		print('\nPathway/RPKM Batch Pathway/Data/Annotation Correlator')
		print('Usage: ')
		print('rpkm_annotate.py <folder containing pathway, data, and annotation files> [output filename] [--select-pathways <file>] [--anno-suffix <suffix>]')
		print('\nIf no output file is specified, defaults to pwy_anno.tsv\n')
		print('The --help flag displays this message.')
		print('Specifying --select-pathways <pathway list file> specifies a file containing a list of pathways to output (ignoring other pathways)')
		print('Specifying --anno-suffix <suffix> sets the filename suffix for the annotation files to be used, allowing for multiple sets of annotation files in the same folder.')
		print('\nSee README for further information.\n')


	target_folder = ""

	select_pathways = False
	pwy_select_filename = ""
	pwys_selected = []
	anno_suffix = '.metacyc-2016-10-31.lastout.parsed.txt'

	args = list(sys.argv)

	# If --help is specified, print usage information and exit.
	if '--help' in args:
		print('Printing usage information.')
		printUsage()
		quit()

	# If --anno-suffix <suffix> is selected, use that as the filename suffix for annotation files
	if '--anno-suffix' in args:
		# Acquire the suffix from the command line arguments
		anno_suffix_idx = args.index('--anno-suffix') + 1
		anno_suffix = args[anno_suffix_idx]
		args.remove(anno_suffix)
		args.remove('--anno-suffix')


	# If --select-pathways <file> is specified, load the list of pathways to process from the specified file
	if '--select-pathways' in args:
		# Acquire the specified file name from the command-line arguments
		select_file_idx = args.index('--select-pathways') + 1
		pwy_select_filename = args[select_file_idx]
		args.remove('--select-pathways')
		args.remove(pwy_select_filename)

		# Read the names of pathways from the specified file
		# Currently, the format of the specified file is hardcoded, and must be as such:
		# - comma-separated
		# - the first row is a header
		# - the first column is the short-names of each pathway to look at
		try:
			with open(pwy_select_filename, 'r') as pwy_sel_file:
				sel_reader = csv.reader(pwy_sel_file, delimiter=',')

				# Keep track of whether the first row is currently being read, so the header can be ignored
				first_row = True

				for row in sel_reader:

					# Ignore the header row
					if first_row == True:
						first_row = False
						continue

					# Append the short-name of the pathway on the current row to the list of pathways to look at
					pwys_selected.append(row[0])

				print('Looking at ' + str(len(pwys_selected)) + ' pathways total.')

		# If an error occurred loading/parsing this file, output the error to the command line and exit. 
		except Exception as e:
			print("Error loading list of selected pathways - exiting.")
			print("Exception: " + str(e))
			quit()


	if len(args) == 2: # Output file not specified
		target_folder = args[1]
		output_filename = 'pwy_anno.tsv'

	elif len(args) == 3: # Output file specified
		target_folder = args[1]
		output_filename = args[2]

	else: # If command-line arguments are not well-formed, print usage information and exit.
		printUsage()
		quit()


	batchCorrelateAnnotate(target_folder, output_filename, selected_pathways=pwys_selected, anno_file_suffix=anno_suffix)
