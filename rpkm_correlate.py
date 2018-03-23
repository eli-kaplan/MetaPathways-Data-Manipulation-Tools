""" rpkm_correlate.py: Library / script to correlate / combine RPKM data values with pathway information """

__author__ = "Eli Kaplan"
__email__ = "eli.kaplan@alumni.ubc.ca"
__license__ = "CC0 - No rights reserved. See LICENSE"

import csv
import sys


def loadPathwayInfoFromFile(filename, csv_separator):
	""" loadPathwayInfoFromFile(): Loads pathway information from given filename (with given CSV separator character)

								   Parameters: 
								   - filename: name of file to load pathway information from (string)
								   - csv_separator: CSV separator character (e.g. '\t')

									Returns: Tuple(String(Sample Name), List(Tuple(Pathway Short Name, Pathway Common Name, List(Pathway ORF IDs))))) """


	# List to place loaded pathway info into
	pathway_info = []

	# Keep track of the name of the sample (i.e. MaxBin_###) to properly match data
	sample_name = ""

	# Load the pathway information from the specified file
	try:
		with open(filename, 'r') as pwy_file:
			pwy_reader = csv.reader(pwy_file, delimiter=csv_separator)

			headerProcessed = False # Keep track of whether the file header has been processed
			
			# Keep track of the sample name to properly match data
			sample_name = ""

			# Column indices for each used field (SAMPLE, PWY_NAME, PWY_COMMON_NAME, ORFS)
			idx_name = -1
			idx_cname = -1
			idx_orfs = -1
			idx_sample = -1

			for row in pwy_reader:
				if headerProcessed == False: # The first line contains a header - use this to find the position of the necessary columns
					# Store the column index for the sample name (SAMPLE), pathway shortened name (PWY_NAME), common name (PWY_COMMON_NAME), and ORF data IDs (ORFS)
					try:
						idx_sample = row.index("SAMPLE")
						idx_name = row.index("PWY_NAME")
						idx_cname = row.index("PWY_COMMON_NAME")
						idx_orfs = row.index("ORFS")
					except:
						# If an index could not be found for all of the required columns, exit. 
						print("ERROR: Pathway file is improperly formatted - not all necessary columns are present.")
						quit()

					headerProcessed = True
					continue
					
				if sample_name == "": # Acquire the sample name if it is not yet set
					sample_name = row[idx_sample]

				pwy_name_str = row[idx_name] 	# PWY_NAME column
				pwy_cname_str = row[idx_cname] 	# PWY_COMMON_NAME column
				pwy_refs_str = row[idx_orfs] 	# ORFS column
				
				# Parse the ORFS column into a list of frame IDs
				pwy_refs_list = pwy_refs_str[1:-1].split(',')

				# Store the parsed (NAME, COMMON_NAME, list(ORFS)) tuple
				pathway_info.append((pwy_name_str, pwy_cname_str, pwy_refs_list))

	except: # If an error occurred while loading and parsing the pathway info file, exit.
		print("ERROR: Could not read/parse pathway information file - exiting.")
		quit()


	# Return the pathway information as a tuple of (Sample Name, list(tuple(Pathway Short Name, Pathway Common Name, list(ORF IDs))
	return (sample_name, pathway_info)




def loadORFDataFromFile(filename, sample_name, csv_separator):
	""" loadORFDataFromFile: Loads ORF RPKM data for a set from a given file (with given sample name and csv separator)

							 Parameters: 
							 - filename: string containing ORF data file
							 - sample_name: string containing name of sample (e.g. MaxBin_33)
							 - csv_separator: string containing CSV separator character (e.g. '\t')

							 Returns: Tuple(Sample Name, List(Tuple(ORF ID, RPKM Reading))) """
	rpkm_data = []

	# Load the RPKM data from the specified file
	try:
		with open(filename, 'r') as data_file:
			data_reader = csv.reader(data_file, delimiter=csv_separator)

			for row in data_reader:
				data_id = 'O_' + row[0].replace(sample_name,'')[1:] # Convert the frame ID column into the O_###_# format used in the pathway info file
				data_value = row[1] 						 		# Capture the data point associated with this frame

				# Store the parsed (ORF ID, Data) tuple
				rpkm_data.append((data_id, data_value))

	except: # If an error occurred while loading / parsing the RPKM data file, exit.
		print("ERROR: Could not read/parse RPKM data file - exiting.")
		quit()

	# Return a tuple of the sample name and the loaded per-ORF RPKM reading data
	return (sample_name, rpkm_data)




def correlatePathwayInfoWithData(sample_name, pathway_info, rpkm_data):
	""" correlatePathwayInfoWtihData: Correlates pathway info and experimental RPKM readings (summing the RPKM readings for each associated ORF) for one sample.

									  Parameters:
									  - sample_name: name of the sample (string)
									  - pathway_info: list of tuples of the following format: (Pathway Short Name, Pathway Common Name, List(ORF IDs as strings))
									  - rpkm_data: list of tuples formatted as: (ORF ID as string, RPKM reading as string (of float))


									  Returns: Tuple(Sample Name, List(Tuple(Pathway Short Name, Pathway Common Name, RPKM Readings Sum)))"""

	# Correlate the data points in pathway_data with the IDs in pathway_info
	pathway_sums = []
	pathway_data_dict = dict(rpkm_data) # Convert the list of id/data tuples to a dict to make searching easier

	# Iterate through the pathways and sum up the ORF data points associated with each ID in the ORFS column
	for pathway in pathway_info:
		pwy_sum = 0.0
		for ref in pathway[2]: # Add the data value associated with each ORF ID to pwy_sum (if it exists)
			if ref in pathway_data_dict:
				pwy_sum += float(pathway_data_dict[ref])
			else:
				print('Missing data point: ' + ref)

		# Store the resulting tuple (PWY_NAME, PWY_COMMON_NAME, sum of referenced data points)
		pathway_sums.append((pathway[0], pathway[1], pwy_sum))

	# Return a tuple of the sample name and the resulting data
	return (sample_name, pathway_sums)




def correlateRPKM(pwy_filename, data_filename, output_filename = 'pwy_data.tsv', csv_separator='\t'):
	""" correlateRPKM(): Takes pathway information (from pwy_filename) and experimental ORF RPKM measurements (from data_filename) and correlates them,
						 summing up RPKM measurements for each pathway by matching ORF IDs between pathways and the RPKM data file. This resulting data 
						 is then output to a file with the specified name (output_filename), which defaults to pwy_data.tsv 

						 Parameters:
						 - pwy_filename: filename of the pathway information file for the sample being analyzed
						 - data_filename: filename of the file containing the RPKM data for each ORF ID
						 - output_filename: file for the results to be stored in
						 - csv_separator: separator in use in the CSV files (e.g. '\t') 

						 Returns: nothing - outputs results to specified file.""" 

	
	# Load the pathway information from the given file
	pwy_info_file_data = loadPathwayInfoFromFile(pwy_filename, csv_separator)
	sample_name = pwy_info_file_data[0] # Name of the sample (E.g. MaxBin_33)
	pathway_info = pwy_info_file_data[1] # List of tuples containing pathway info (short name, common name, list of ORF IDs)

	print("Loaded sample data for " + sample_name + " from " + pwy_filename)


	# Load the ORF RPKM data from the given file
	pwy_orf_file_data = loadORFDataFromFile(data_filename, sample_name, csv_separator)
	pathway_data = pwy_orf_file_data[1]

	print("Loaded RPKM data for " + sample_name + " from " + data_filename)


	# Correlate the pathway information with RPKM data
	pwy_correlation_data = correlatePathwayInfoWithData(sample_name, pathway_info, pathway_data)
	pathway_sums = pwy_correlation_data[1] # Tuple(PWY_NAME, PWY_COMMON_NAME, sum of referenced RPKM data points for this sample)


	# Output the resulting data to the specified file
	try:
		with open(output_filename, 'w') as output_file:
			output_writer = csv.writer(output_file, delimiter=csv_separator)

			for pair in pathway_sums:
				output_writer.writerow(pair)
	
	except: # If an error occurred while writing out the results, exit.
		print("ERORR: File output failed - exiting.")
		quit()

	print("Correlated data for " + sample_name + " output to " + output_filename)




# Allow the script to be run stand-alone (and prevent the following code from running when imported)
if __name__ == "__main__":

	def printUsage():
		"""Prints usage information for this script."""
		print('\nPathway/RPKM Data Correlator')
		print('Usage: ')
		print('rpkm_correlate.py <pathway info file> <RPKM data file> [output file]')
		print('If no output file is specified, defaults to pwy_data.tsv')
		print('\nSee README for further information.\n')

	# Acquire input / output filenames from command-line arguments, if they are well-formed
	info_filename = ""
	data_filename = ""
	output_filename = ""

	if len(sys.argv) == 3: # Output file not specified
		info_filename = sys.argv[1]
		data_filename = sys.argv[2]
		output_filename = 'pwy_data.tsv'

	elif len(sys.argv) == 4: # Output file specified
		info_filename = sys.argv[1]
		data_filename = sys.argv[2]
		output_filename = sys.argv[3]

	else: # If command-line arguments are not well-formed, print usage information and exit.
		printUsage()
		quit()


	correlateRPKM(info_filename, data_filename, output_filename)
