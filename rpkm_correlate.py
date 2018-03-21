""" rpkm_correlate.py: Script to correlate / combine RPKM data values with pathway information """

__author__ = "Eli Kaplan"
__email__ = "eli.kaplan@alumni.ubc.ca"
__license__ = "CC0 - No rights reserved. See LICENSE"

import csv
import sys

def correlateRPKM(pwy_filename, data_filename, output_filename = 'pwy_data.tsv', csv_separator='\t'):
	""" correlateRPKM(): Takes pathway information (from pwy_filename) and experimental ORF RPKM measurements (from data_filename) and correlates them,
						 summing up RPKM measurements for each pathway by matching ORF IDs between pathways and the RPKM data file. This resulting data 
						 is then output to a file with the specified name (output_filename), which defaults to pwy_data.tsv """ 

	pathway_info = []
	pathway_data = []

	# Keep track of the name of the sample (i.e. MaxBin_###) to properly match data
	sample_name = ""


	# Load the pathway information from the specified file
	try:
		with open(pwy_filename, 'r') as pwy_file:
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

	print("Loaded sample data for " + sample_name + " from " + pwy_filename)


	# Load the RPKM data from the specified file
	try:
		with open(data_filename, 'r') as data_file:
			data_reader = csv.reader(data_file, delimiter=csv_separator)

			for row in data_reader:
				data_id = 'O_' + row[0].replace(sample_name,'')[1:] # Convert the frame ID column into the O_###_# format used in the pathway info file
				data_value = row[1] 						 		# Capture the data point associated with this frame

				# Store the parsed (ORF ID, Data) tuple
				pathway_data.append((data_id, data_value))

	except: # If an error occurred while loading / parsing the RPKM data file, exit.
		print("ERROR: Could not read/parse RPKM data file - exiting.")
		quit()


	print("Loaded RPKM data for " + sample_name + " from " + data_filename)


	# Correlate the data points in pathway_data with the IDs in pathway_info
	pathway_sums = []
	pathway_data_dict = dict(pathway_data) # Convert the list of id/data tuples to a dict to make searching easier

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


	# Output the resulting data to the specified file
	try:
		with open(output_filename, 'wb') as output_file:
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
		print('If no output file is specified, defaults to pwy_data.tsv\n')


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