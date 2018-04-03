""" rpkm_correlate_batch.py: Library / script to correlate / combine RPKM data values from multiple sample sets with pathway information """

__author__ = "Eli Kaplan"
__email__ = "eli.kaplan@alumni.ubc.ca"
__license__ = "CC0 - No rights reserved. See LICENSE"

import csv
import sys

from os import listdir
from os.path import join

from rpkm_correlate import loadPathwayInfoFromFile, loadORFDataFromFile, correlatePathwayInfoWithData

def batchCorrelateRPKM(file_dir, output_filename = 'pwy_data_batch.tsv', csv_separator='\t', pwy_file_suffix='.pwy.txt', data_file_suffix='.orf_rpkm.txt', excl_zeroes=False, stats_file_suffix='_stats', separate_stats = True):
	""" batchCorrelateRPKM(): Takes a directory with two sets of files (one set of files containing pathway information and assocated ORF RPKM data IDs, and the
							  other containing ORF IDs and associated data points), sums the RPKM data for each pathway within a sample, and combines all of the sample results
							  into a single file.

							  Parameters:
							  - file_dir: path to directory containing necessary files
							  - output_filename: name of file to output results to 
							  - csv_separator: column separator used in given files
							  - pwy_file_suffix: suffix for files containing pathway information
							  - data_file_suffix: suffix for files containing RPKM data 
							  - excl_zeroes: whether average calculations should exclude zero-values (i.e. pathways not found in specific samples)

							  Returns: nothing - outputs results to specified file"""

	# Create a list of all files in the target directory
	all_files = [join(file_dir, f) for f in listdir(file_dir)]

	# Split the list of all the files into two lists: a list of pathway information files, and a list of RPKM data files
	pwy_files = []
	data_files = []

	for file in all_files:
		if pwy_file_suffix in file:
			pwy_files.append(file)
		elif data_file_suffix in file:
			data_files.append(file)
		else:
			print("Unknown file in batch directory: " + file)


	# Match each pathway information file with its corresponding RPKM data file (assuming they have the same sample name)
	file_pairs = []

	for pwy_file in pwy_files:
		# Generate the name of the corresponding data file based on the name of the pathway info file
		corresponding_data_file = pwy_file.replace(pwy_file_suffix, data_file_suffix)
	
		# If this file exists, pair it with the pathway file
		if corresponding_data_file in data_files:
			file_pairs.append((pwy_file, corresponding_data_file))

		# Otherwise, print a warning and don't add the anything to the list of pairs of files to work with.
		else:
			print("Missing ORF data file: " + corresponding_data_file)


	# Load/parse each pathway information file and corresponding data file, then correlate each pathway with its respective data points
	all_files_data = []

	for (pathway_file, data_file) in file_pairs:
		# Load pathway information from the given file
		pathway_info = loadPathwayInfoFromFile(pathway_file, csv_separator)
		cur_sample = pathway_info[0]

		# Load RPKM data from the file corresponding to the pathway info file
		rpkm_data = loadORFDataFromFile(data_file, cur_sample, csv_separator)

		# Correlate the data from the two files
		corr_pwy_data = correlatePathwayInfoWithData(cur_sample, pathway_info[1], rpkm_data[1])

		# Store this correlated data (tuple of (Sample Name, Data)) in all_files_data
		all_files_data.append(corr_pwy_data)

		print("Loaded data for sample: " + cur_sample)				


	# Create a dict containing all of this previously loaded data, keyed by pathway short-name (such as to allow per-sample data to be accessed for each pathway)
	per_pathway_data = dict()
	all_samples = [] # Keep track of all of the samples encountered

	# Add each pathway from each sample to per_pathway_data
	for (sample, pathways) in all_files_data:
		for pwy in pathways:
			pwy_name = pwy[0] # Short name for this pathway
			pwy_cname = pwy[1] # Common name for this pathway
			pwy_rpkm = pwy[2] # Sum of RPKM readings for this pathway for this sample

			# If the pathway is not in the dict yet, add it.
			if pwy_name not in per_pathway_data:
				per_pathway_data[pwy_name] = [pwy_name, pwy_cname, dict()] # the dict() field will contain {'Sample Name' : 'RPKM Sum'} data for each pathway
			
			# Append the data for this current sample to the data for this pathway
			per_pathway_data[pwy_name][2][sample] = pwy_rpkm

		# Add this sample to the list of samples encountered
		all_samples.append(sample)


	# For each pathway, if said pathway is not measured in one of the samples, assign a RPKM sum value of 0.0 for that sample 
	for pathway, data in per_pathway_data.items():
		for sample in all_samples:
			if sample not in data[2]:
				data[2][sample] = 0.0

	
	# Generate a header for the output tabulated file 
	output_file_header = ['Name', 'Common Name', 'Average RPKM', 'RPKM Sum', 'In # Samples']
	for sample in sorted(all_samples):
		output_file_header.append(sample)


	# Keep track of the per-sample RPKM sums and numbers of non-zero values as {'Sample Name' : #}
	sample_col_sums = {}
	sample_col_nonzero_values = {}

	# Output the resulting data to the chosen file
	try:
		with open(output_filename, 'w') as output_file:
			output_writer = csv.writer(output_file, delimiter=csv_separator)

			# Write the header to the output file as the first line.
			output_writer.writerow(output_file_header)


			# Write out the data for each pathway
			for pathway, data in per_pathway_data.items():
				# Add the pathway short-name and common-name to the current row to be written
				row = [data[0], data[1]]

				# Calculate the per-pathway RPKM value sum across all of the samples
				rpkm_sum = sum(data[2].values())


				# Calculate the per-pathway RPKM value average
				rpkm_average = 0

				# Only perform this calculation for pathways that have been measured in at least one sample
				if rpkm_sum > 0.0:

					# If zero-values are to be excluded, calculate the number of samples this pathway was measured in
					if excl_zeroes == True:
						num_nonzero_samples = 0
						for sample_reading in data[2].values():
							if sample_reading != 0.0:
								num_nonzero_samples += 1
						
						# Calculate the average as (sum) / (number of samples where this pathway has been identified)
						if num_nonzero_samples != 0:
							rpkm_average = rpkm_sum / num_nonzero_samples

					# Otherwise, just use the total number of samples to calculate the average
					else:
						rpkm_average = rpkm_sum / len(data[2].values())

				# Append the per-pathway RPKM average and sum to the row
				row.append(rpkm_average)
				row.append(rpkm_sum)

				# Calculate fraction of samples that the pathway appears in
				in_n_samples = 0
				for sample, val in data[2].items():
					if val != 0.0:
						in_n_samples += 1

				# Append this to the current row to be written 
				row.append(str(in_n_samples) + '/' + str(len(all_samples)))

				# Iterate over all of the {'Sample' : 'RPKM Sum'} data for this pathway
				for sample, val in sorted(data[2].items()):

					# Add the reading for this pathway in each sample to the dict of per-sample total RPKM sums
					if sample in sample_col_sums.keys():
						sample_col_sums[sample] += val
					else:
						sample_col_sums[sample] = val

					# If the reading for this pathway is non-zero, add it to the dict of per-sample number of non-zero-RPKM pathways (i.e. unique pathways found in each sample)
					if val != 0.0:
						if sample in sample_col_nonzero_values.keys():
							sample_col_nonzero_values[sample] += 1
						else:
							sample_col_nonzero_values[sample] = 1


					# Append the RPKM reading for this pathway in this sample to the row to be written out 
					row.append(val)

				# Write the current row to the output file
				output_writer.writerow(row)


			# If per-sample stats are not specified to be separated, place them in the bottom of the output file
			if separate_stats == False:
				# Add a row for the total per-sample RPKM sums to the bottom of the file
				sums_row = ['SAMPLE-SUMS', 'Per-Sample RPKM Sum', '--', '--', '--']	

				for sample, total_sum in sorted(sample_col_sums.items()):
					sums_row.append(total_sum)

				output_writer.writerow(sums_row)


				# Add a row for the total per-sample RPKM averages to the bottom of the file 
				averages_row = ['SAMPLE-AVGS', 'Per-Sample RPKM Average', '--', '--', '--']

				total_num_pathways = len(per_pathway_data.keys())

				# If zeroes are excluded, calculate this average as (per-sample RPKM sum) / (per-sample number of unique pathways observed with non-zero RPKM)
				if excl_zeroes == True:
					for sample, total_sum in sorted(sample_col_sums.items()):
						averages_row.append(total_sum / sample_col_nonzero_values[sample])

				# Otherwise, calculate as (per-sample RPKM sum) / (total number of pathways loaded from all files)
				else:
					for sample, total_sum in sorted(sample_col_sums.items()):
						averages_row.append(total_sum / total_num_pathways)

				output_writer.writerow(averages_row)


		# If the per-sample stats are supposed to be separated, create a new file for them and output them there.
		if separate_stats == True:
			with open(output_filename + stats_file_suffix, 'w') as output_stats_file:
				stats_writer = csv.writer(output_stats_file, delimiter=csv_separator)

				# Write out the header for this file
				stats_header = ['Sample', 'Total Sample RPKM Sum', 'Average Per-Pathway RPKM']
				stats_writer.writerow(stats_header)

				# Write out a statistics row for each sample
				for sample, total_sum in sample_col_sums.items():
					cur_row = [sample, total_sum]

					if excl_zeroes == True:
						cur_row.append(total_sum / sample_col_nonzero_values[sample])

					else:
						total_num_pathways = len(per_pathway_data.keys())
						cur_row.append(total_sum / total_num_pathways)

					stats_writer.writerow(cur_row)


	
	except Exception as e: # If an error occurred while writing out the results, exit.
		print("ERORR: File output failed - exiting.")
		print("Exception: " + str(e))
		quit()


	




# Allow the script to be run stand-alone (and prevent the following code from running when imported)
if __name__ == "__main__":

	def printUsage():
		"""Prints usage information for this script."""
		print('\nPathway/RPKM Batch Data Correlator')
		print('Usage: ')
		print('rpkm_correlate_batch.py <pathway info / RPKM data files folder> [output file] [--exclude-zeroes] [--separate-stats]')
		print('\nIf no output file is specified, defaults to pwy_data_batch.tsv\n')
		print('The --exclude-zeroes flag excludes zero-values from all average calculations.')
		print('The --help flag shows this information.')
		print('\nSee README for further information.\n')


	target_folder = ""

	args = list(sys.argv)
	exclude_zeroes = False
	sep_stats = False

	# If --exclude-zeroes is specified, enable that flag in batchCorrelateRPKM()
	if '--exclude-zeroes' in args:
		exclude_zeroes = True
		args.remove('--exclude-zeroes')
		print('Excluding zero-values from average calculation.')


	# If --separate-stats is specified, output per-sample stats to another file
	if '--separate-stats' in args:
		sep_stats = True
		args.remove('--separate-stats')
		print('Will output per-sample stats to a separate file.')


	# If --help is specified, print usage information and exit.
	if '--help' in args:
		print('Printing usage information.')
		printUsage()
		quit()


	if len(args) == 2: # Output file not specified
		target_folder = args[1]
		output_filename = 'pwy_data_batch.tsv'

	elif len(args) == 3: # Output file specified
		target_folder = args[1]
		output_filename = args[2]

	else: # If command-line arguments are not well-formed, print usage information and exit.
		printUsage()
		quit()


	batchCorrelateRPKM(target_folder, output_filename, excl_zeroes=exclude_zeroes, separate_stats=sep_stats)
