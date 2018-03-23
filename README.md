# Pathway RPKM Correlation Script / Library
- Takes metabolic pathway information and Open Reading Frame (ORF) Reads Per Kilobase Million (RPKM) data from separate files, correlates RPKM data points with their respective pathways, and manipulates / outputs this data.


# Usage

## rpkm_correlate

### As a Script
```
python rpkm_correlate.py (pathway information file) (RPKM data file) [output file]
```
If not specified, the output file defaults to ``` pwy_data.tsv ```, in the working directory. 

Running this script will correlate the RPKM data points from the given data file with their respective pathways in the given pathway information file (based on ORF ID), sum the RPKM readings for each pathway's frames in this sample, and output this combined data to the specified output file (or ```pwy_data.tsv```). Files must be formatted as specified in ```File Format``` for this script to function properly.

### As a Python Library
```
from rpkm_correlate.py import correlateRPKM

correlateRPKM("pathway_info.txt", "rpkm_data.txt", "pwy_data.tsv")
```
This does exactly what running the script (```As a Script``` above) would do. The library also includes the following other functions, useful for working with these types of files (see code for documentation):
- ```loadPathwayInfoFromFile(filename, csv_separator)```
- ```loadORFDataFromFile(filename, sample_name, csv_separator)```
- ```correlatePathwayInfoWithData(sample_name, pathway_info, rpkm_data)```

The ```correlateRPKM()``` function also has an optional parameter, ```csv_separator```, which specifies the separator character used in the input and output files (defaults to ```'\t'```, a tab character). If the input files are in CSV (comma-separated) format, as opposed to TSV (tab-separated), add the following parameter: 
```csv_separator = ','```. 


## rpkm_correlate_batch
### As a Script
``` 
python rpkm_correlate_batch.py <pathway info / RPKM data files folder> [output file] [--exclude-zeroes] [--separate-stats]
```
Running this script takes a directory full of pathway information and RPKM data point files (see ```File Format```), and calculates the sum of all of the RPKM readings for each pathway within each sample. It then produces a new file (```pwy_data_batch.tsv``` if not specified), which contains a tabulated listing of every pathway (in all of the files in the specified directory) and the per-sample RPKM reading sum for each of these pathways. Additionally, the sum of all of the RPKM readings for each pathway and for each sample, as well as the per-sample and per-pathway RPKM averages, are also stored in the file. The per-pathway calculations are stored in two columns, and the per-sample calculations are stored on the last two rows of the file. 

Three command-line options can be specified to this script:
- ```--exclude-zeroes``` excludes zero-values from the calculation of averages. This allows averages to be calculated for each pathway only for samples it is present in, and allows averages for each sample to be calculated to only include pathways actually present in that sample.
- ```--separate-stats``` instructs the program to output per-sample stats in a separate file (instead of on the last two lines of the output file). This file will be named the same as the output file, but with ```_stats``` appended to the name.
- ```--help``` prints usage information to the console.


### As a Python Library
```
from rpkm_correlate_batch import batchCorrelateRPKM

batchCorrelateRPKM("input_directory/", output_filename = 'pwy_data_batch.tsv', separate_stats=False)
```
This will perform the same actions as the default usage of the script. There are a few optional parameters that can be passed to this function:
- ```csv_separator``` changes the separator character used to tabulate the input/output files (default `\t`)
- ```pwy_file_suffix``` is the suffix for pathway information files (to figure out which files they are in the specified directory). Defaults to `.pwy.txt`
- ``` data_file_suffix``` is the suffix for RPKM data files. Defaults to ```.orf_rpkm.txt```
- ``` excl_zeroes``` determines whether zero-values are excluded from the average calculations.
- ```separate_stats``` determines whether per-sample statistics are placed in the main output file, or a separate file. Defaults to `True` (separate file).

## File Format 
The input files should be formatted as TSV/CSV files, with the following requirements:
- The pathway information file should have the suffix `.pwy.txt`, have at least 4 columns, with a header on the first line: 
	- ```SAMPLE``` (name of the sample, e.g. ```MaxBin27```)
	- ```PWY_NAME``` (short-name for the pathway, e.g. ```PWY-4416```)
	- ```PWY_COMMON_NAME``` (common-name for the pathway, e.g. ```CMP phosphorylation```)
	- ```ORFS``` (list of frame IDs, e.g. ```[O_7_7,O_164_9]```)


- The RPKM data files should have the prefix `.orf_rpkm.txt`, and should have two columns, with no header:
	- The first column should be the ORF ID, formatted as ```[Sample Name]_###_##```
		- For example, if the sample name is ```MaxBin_22```, and the frame ID is ```O_145_2```, the resulting ID would be ```MaxBin22_145_2```
	- The second column should be the RPKM data point, formatted as a ```float``` (e.g. ```1.9415```)

For `rpkm_correlate_batch` to function correctly, the name for each file in each pathway info / data file pair must match except for the prefix. If the pathway information file is named `maxbin_44.pwy.txt`, the data file must be named `maxbin_44.orf_rpkm.txt` - this is case-sensitive!


# License
This work is available under the ```Creative Commons CC0 1.0 Universal license```. You can modify, distribute, and otherwise use this work, even for commercial purposes, without permission or notice. See LICENSE for full details. 

If you find this tool helpful, feel free to let me know! (Although, per the license, you aren't required to.) 



