# Pathway RPKM Correlation Script / Library
- Takes metabolic pathway information and Open Reading Frame (ORF) Reads Per Kilobase Million (RPKM) data from separate files, correlates RPKM data points with each pathway, sums the RPKM values per pathway, and outputs to a new file. 



## Usage

### As a Script
```
python rpkm_correlate.py (pathway information file) (RPKM data file) [output file]
```
If not specified, the output file defaults to ``` pwy_data.tsv ```, in the working directory. 


### As a Python Library
```
from rpkm_correlate.py import correlateRPKM

correlateRPKM("pathway_info.txt", "rpkm_data.txt", "output.txt")
```
The ```correlateRPKM()``` function also has an optional parameter, ```csv_separator```, which specifies the separator character used in the input and output files (defaults to ```'\t'```, a tab character). If the input files are in CSV (comma-separated) format, as opposed to TSV (tab-separated), add the following parameter: ```csv_separator = ','```. 


### File Format 
The input files should be formatted as TSV/CSV files (depending on the separator), with the following requirements:
- The pathway information file should have at least 4 columns, with a header on the first line: 
	- ```SAMPLE``` (name of the sample, e.g. ```MaxBin27```)
	- ```PWY_NAME``` (short-name for the pathway, e.g. ```PWY-4416```)
	- ```PWY_COMMON_NAME``` (common-name for the pathway, e.g. ```CMP phosphorylation```)
	- ```ORFS``` (list of frame IDs, e.g. ```[O_7_7,O_164_9]```)


- The RPKM data file should have two columns, with no header:
	- The first column should be the ORF ID, formatted as ```[Sample Name]_###_##```
		- For example, if the sample name is ```MaxBin_22```, and the frame ID is ```O_145_2```, the resulting ID would be ```MaxBin22_145_2```
	- The second column should be the RPKM data point, formatted as a ```float``` (e.g. ```1.9415```)


## License
This work is available under the ```Creative Commons CC0 1.0 Universal license```. You can modify, distribute, and otherwise use this work, even for commercial purposes, without permission or notice. See LICENSE for full details. 

