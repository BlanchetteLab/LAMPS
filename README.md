# LAMPS
Sequence analysis pipeline for 2C-ChIP and 5C products

## Overview
The 'Ligation-mediated Amplified, Multiplexed Paired-end Sequence' or LAMPS is is a Linux/MacOS command line interface for analyzing paired-end sequences, which may or may not be multiplexed.

## Software requirements
1) Python (v2.7.15 or v3.8.0 tested): https://conda.io/docs/user-guide/install/download.html (recommended)
2) BLAST (v2.5.0+ tested): https://www.ncbi.nlm.nih.gov/books/NBK279671/
3) SAMtools (v1.3.1 tested - optional for BAM file processing): via Anaconda (https://conda.io/docs/user-guide/install/download.html) or https://formulae.brew.sh/formula/samtools or http://samtools.sourceforge.net/
4) *MacOS users only* - gnu-sed (v4.7 tested): https://formulae.brew.sh/formula/gnu-sed

## Installation guide
Installation is expected to take a few minutes:
1) Either download the package by clicking the "Clone or download" button, unzipping file in desired location, and renaming the directory "LAMPS"   OR   Use the command line ``` git clone https://github.com/BlanchetteLab/LAMPS ```.
2) If any of the required Python modules are not installed, install them using Anaconda (https://conda.io/docs/user-guide/install/download.html).
3) If BLAST is not installed, install it (https://www.ncbi.nlm.nih.gov/books/NBK279671/) and make sure the BLAST executable is in your path.
5) If sequencing files are in BAM format and SAMtools is not installed, install it (http://www.htslib.org/download/) and make sure the SAMtools executable is in your path.
6) Download and uncompress the example directories for the following test 2C-ChIP and 5C data sets from Wang et al. 2019 into the LAMPS/ directory:
    * 2C-ChIP example: https://www.cs.mcgill.ca/~blanchem/LAMPS/2C-ChIP.tar.gz
    * 5C example: https://www.cs.mcgill.ca/~blanchem/LAMPS/5C.tar.gz

## Quick start

Set current directory to LAMPS:

```cd LAMPS/```

Process 2C-ChIP sequencing data:

```python LAMPS.py ./2C-ChIP/LAMPS.config ./2C-ChIP/HOXA_2C-ChIP_primers.txt 2C-ChIP ./2C-ChIP/LAMPS_output```

Process 5C sequencing data:

```python LAMPS.py ./5C/LAMPS.config ./5C/HOXA_5C_primers.txt 5C ./5C/LAMPS_output```

## Input
LAMPS config file - human-readable text file (tab-separated values [TSV] format) with three required and 5 optional columns:
1) *required*: paired-end library name
2) *required*: barcode sequence for multiplexing - if no barcode, provide an empty string
3) *required*: complete file path to sequencing file (FASTQ or BAM)
4) *optional*: primers to exclude - comma separated list of primer names found in the second column of the primer file (**must be identical**)
5) *optional*: batch number - user-assigned batch id that groups libraries with respective input library for normalization
6) *optional*: normalization factors - one or more normalization factors found in their own additional columns (i.e., if there are three normalization factors, then there will be 8 columns present in the config file). LAMPS will multiply factors together.

&nbsp;&nbsp;&nbsp;&nbsp;*Columns 5-8 are used when normalizing 2C-ChIP libraries and not typically applied to 5C data

&nbsp;&nbsp;&nbsp;&nbsp;*Example config files can be found in the example 2C-ChIP and 5C folders

Primer file - human-readable text file (TSV format) with eight required columns:
1) *required*: genome locus name
2) *required*: primer name - for excluding primer pairs, these column values must be identical to column four's in the config file
3) *required*: primer strand - forward 'F' or reverse 'R'
4) *required*: primer sequence (**without barcode**)
5) *required*: genome assembly
6) *required*: chromosome
7) *required*: primer starting basepair
8) *required*: primer ending basepair

&nbsp;&nbsp;&nbsp;&nbsp;*Example primer files can be found in the example 2C-ChIP and 5C folders

## Output

results/ - directory containing the following final output of LAMPS:
* *.raw.matrix - contains the frequency count of target sequences found within the sequenced library. Note - rows/cols are FWD then RVS primers resulting in the following four quadrants: F-F,F-R,R-F,R-R
* *log_raw.heatmap.png - heat map of the log(values) found in *.raw.matrix

&nbsp;&nbsp;&nbsp;&nbsp;*2C-ChIP specific output in the results/*:
* *.raw.bedGraph - raw frequency count of on-diagonal primer pairs in bedGraph format (easily viewable on most genome browsers)
* *.norm.bedGraph - normalized frequency count of on-diagonal primer pairs (will be input normalized if input library was found - please see LAMPS stdout to be sure, a warning will be thrown if not found)
* *.line_plot.png - line plots of either raw or normalized frequency counts of on-diagonal primer pairs
* 'raw_totals_vs_norm_factor.scatter.png' - scatter plot of normalization factor vs. raw total counts of libraries (expected to show a linear trend - Spearman correlation provided)
* 'raw_totals_vs_norm_factor.tsv' - values of normalization factor vs. raw total counts of libraries

### Secondary output

The following directories contain files specific to primer or mapping Quality Control (QC).

primer_files/ - directory containing the following files relevant to the primer pair analysis:
* *.fasta - ligated primer pair and short read FASTA files used as input to BLAST
* \*.matrix and *.png - primer pair similarities (bitscore between sequences) matrix and heatmap. Useful for identifying problematic primers.

mapping_files/ - directory containing the following BLASTn mapping and summary files:
* BLAST/ - directory containing custom BLASTdb files
* \*.fasta (*.fastq if conversion from BAM required) - FASTA formatted files of sequencing data used as input to BLAST
* *.BLAST.tsv - BLASTn output of mapping read data against the ligated primer pair library
* *.BLAST.log - BLASTn stdout from mapping
* *.BLAST_summary.bar_plot.png - plot of total and mapped reads to the ligated primer pair library
* 'BLAST_summary.word_size_<>.tsv' - counts of total and mapped reads for libraries
* *.read_count.bar_plot.png -  bar plot breakdown of sequenced reads
* 'read_count_frequency_table.tsv' - counts of sequenced reads breakdown
* short_read_analysis/ - contains the following files relevant to the secondary mapping of 'unmappables' to the short read library:
     * *.BLAST.tsv - BLASTn output of mapping read data against the short read library
     * *.BLAST.log - BLASTn stdout from mapping
     * *.BLAST_summary.bar_plot.png - plot of total and mapped reads to the short read library
     * *.short_reads_summary.tsv - counts of total and mapped reads for library

## Command line details:
Processing sequenced 2C-ChIP or 5C data:
```
LAMPS.py [-h] [--num_cpus NUM_CPUS] config primers type output

positional arguments:
  config			file path to LAMPS config file
  primers			file path to TSV file containing primer information
  type				source of paired-end reads:[2C-ChIP,5C]
  output			file path to output folder	

optional arguments:
  -h, --help		show this help message and exit
  --num_cpus		set the number of cpus - default = total_num_cpus-2
  --word_size     set the minimum required sequence length for processing - default = automated calculation
```

## Testing
All software was tested on Linux Ubuntu 12.04.5 LTS (GNU/Linux 3.2.0-86-generic x86_64).

## License:
LAMPS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation.

LAMPS is distributed in the hopes that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
