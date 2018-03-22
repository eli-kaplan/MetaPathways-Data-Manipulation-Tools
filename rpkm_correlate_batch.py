""" rpkm_correlate_batch.py: Library / script to correlate / combine RPKM data values from multiple sample sets with pathway information """

__author__ = "Eli Kaplan"
__email__ = "eli.kaplan@alumni.ubc.ca"
__license__ = "CC0 - No rights reserved. See LICENSE"

import csv
import sys

from os import listdir
from os.path import join

from rpkm_correlate import loadPathwayInfoFromFile, loadORFDataFromFile, correlatePathwayInfoWithData

def batchCorrelateRPKM(file_dir, output_filename = 'pwy_data_batch.tsv', csv_separator='\t', pwy_file_suffix='.pwy.txt', data_file_suffix='.orf_rpkm.txt'):
	""" batchCorrelateRPKM(): Takes a directory with two sets of files (one set of files containing pathway information and assocated ORF RPKM data IDs, and the
							  other containing ORF IDs and associated data points), sums the RPKM data for each pathway within a sample, and combines all of the sample results
							  into a single file.

							  Parameters:
							  - file_dir: path to directory containing necessary files
							  - output_filename: name of file to output results to 
							  - csv_separator: column separator used in given files
							  - pwy_file_suffix: suffix for files containing pathway information
							  - data_file_suffix: suffix for files containing RPKM data 

							  Returns: nothing - outputs results to specified file"""

	all_files = [join(file_dir, f) for f in listdir(file_dir)]

	pwy_files = []
	data_files = []

	for file in all_files:
		if pwy_file_suffix in file:
			pwy_files.append(file)
		elif data_file_suffix in file:
			data_files.append(file)
		else:
			print("Unknown file in batch directory: " + file)


	file_pairs = []

	for pwy_file in pwy_files:
		corresponding_data_file = pwy_file.replace(pwy_file_suffix, data_file_suffix)
	
		if corresponding_data_file in data_files:
			file_pairs.append((pwy_file, corresponding_data_file))

		else:
			print("Missing ORF data file: " + corresponding_data_file)


	all_files_data = []

	for (pathway_file, data_file) in file_pairs:
		pathway_info = loadPathwayInfoFromFile(pathway_file, csv_separator)
		cur_sample = pathway_info[0]
		rpkm_data = loadORFDataFromFile(data_file, cur_sample, csv_separator)

		corr_pwy_data = correlatePathwayInfoWithData(cur_sample, pathway_info[1], rpkm_data[1])

		all_files_data.append(corr_pwy_data)

		print("Loaded data for sample: " + cur_sample)				


	per_pathway_data = dict()
	all_samples = []

	for (sample, pathways) in all_files_data:
		for pwy in pathways:
			pwy_name = pwy[0] # Short name for this pathway
			pwy_cname = pwy[1] # Common name for this pathway
			pwy_rpkm = pwy[2] # Sum of rpkm readings for this pathway for this sample

			if pwy_name not in per_pathway_data:
				per_pathway_data[pwy_name] = [pwy_name, pwy_cname, dict()]
			
			per_pathway_data[pwy_name][2][sample] = pwy_rpkm

		all_samples.append(sample)


	# For each pathway, if said pathway is not measured in one of the samples, assign a RPKM sum value of 0.0 for that sample 
	for pathway, data in per_pathway_data.iteritems():
		for sample in all_samples:
			if sample not in data[2]:
				data[2][sample] = 0.0


	#for pathway, data in per_pathway_data.iteritems():
	#	print(pathway + " : " + str(data[2]))
	

	all_samples_dict = {}
	for sample in all_samples:
		all_samples_dict[sample] = ""

	all_samples_sorted = []
	for sample, unused in all_samples_dict.iteritems():
		all_samples_sorted.append(sample)

	
	# Generate a header for the output CSV file 
	output_file_header = ['Name', 'Common Name', 'Average RPKM', 'RPKM Sum']
	for sample in all_samples_sorted:
		output_file_header.append(sample)




	# Output the resulting data to the chosen file
	try:
		with open(output_filename, 'w') as output_file:
			output_writer = csv.writer(output_file, delimiter=csv_separator)

			# Write the header to the output file
			output_writer.writerow(output_file_header)

			# Write out the pathway data
			for pathway, data in per_pathway_data.iteritems():
				row = [data[0], data[1]]

				rpkm_sum = sum(data[2].values())

				# TODO: flag to exclude zeroes from RPKM averages
				rpkm_average = 0

				if rpkm_sum > 0:
					rpkm_average = rpkm_sum / len(data[2].values())

				row.append(rpkm_average)
				row.append(rpkm_sum)

				for sample, val in data[2].iteritems():
					row.append(val)

				output_writer.writerow(row)

	
	except: # If an error occurred while writing out the results, exit.
		print("ERORR: File output failed - exiting.")
		quit()


	




# Allow the script to be run stand-alone (and prevent the following code from running when imported)
if __name__ == "__main__":

	def printUsage():
		"""Prints usage information for this script."""
		print('\nPathway/RPKM Batch Data Correlator')
		print('Usage: ')
		print('rpkm_correlate_batch.py <pathway info / RPKM data files folder> [output file]')
		print('If no output file is specified, defaults to pwy_data_batch.tsv\n')


	target_folder = ""

	if len(sys.argv) == 2: # Output file not specified
		target_folder = sys.argv[1]
		output_filename = 'pwy_data_batch.tsv'

	elif len(sys.argv) == 4: # Output file specified
		target_folder = sys.argv[1]
		output_filename = sys.argv[2]

	else: # If command-line arguments are not well-formed, print usage information and exit.
		printUsage()
		quit()


	batchCorrelateRPKM(target_folder, output_filename)
