# LAMPS
Sequence analysis pipeline for 2C-ChIP and 5C products

## Overview
The 'Ligation-mediated Amplified, Multiplexed Primer-pair Sequence' or LAMPS is is a Linux/MacOS command line interface for analyzing Ligation-Mediated Amplification (LMA) sequences, which may or may not be multiplexed.

## Software requirements
1) Python (v2.7.15 or v3.8.1 tested): https://conda.io/docs/user-guide/install/download.html (recommended)
2) Bowtie 2 (v2.3.4.2 tested): http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2 (recommended)\
OR\
BLAST (v2.5.0+ tested): https://www.ncbi.nlm.nih.gov/books/NBK279671/
4) SAMtools (v1.3.1 tested - optional for BAM file processing): via Anaconda (https://conda.io/docs/user-guide/install/download.html) or https://formulae.brew.sh/formula/samtools or http://samtools.sourceforge.net/
5) *MacOS users only* - gnu-sed (v4.7 tested): https://formulae.brew.sh/formula/gnu-sed

## Installation guide
Installation is expected to take a few minutes:
1) Either download the package by clicking the "Clone or download" button, unzipping file in desired location, and renaming the directory "LAMPS"   OR   Use the command line ``` git clone https://github.com/BlanchetteLab/LAMPS ```.
2) If any of the required Python modules are not installed, install them using Anaconda (https://conda.io/docs/user-guide/install/download.html).
3) If Bowtie 2 or BLAST are not installed, install one and make sure the executable is in your path.
5) If sequencing files are in BAM format and SAMtools is not installed, install it (http://www.htslib.org/download/) and make sure the SAMtools executable is in your path.
6) Download and uncompress the example directories for the following test 2C-ChIP and 5C data sets from Wang et al. 2019 into the LAMPS/ directory:
    * 2C-ChIP example: https://www.cs.mcgill.ca/~blanchem/LAMPS/2C-ChIP.tar.gz
    * 5C example: https://www.cs.mcgill.ca/~blanchem/LAMPS/5C.tar.gz

## Quick start

Set current directory to LAMPS:

```cd LAMPS/```

Process 2C-ChIP sequencing data:

```python LAMPS.py ./2C-ChIP/LAMPS.config ./2C-ChIP/HOXA_2C-ChIP_primers.txt 2C-ChIP ./2C-ChIP/LAMPS_output```

Note - creation of Bowtie 2 indices for 2C-ChIP ligated primer sequences takes ~1.5 hours to complete

Process 5C sequencing data:

```python LAMPS.py ./5C/LAMPS.config ./5C/HOXA_5C_primers.txt 5C ./5C/LAMPS_output```

## Input
LAMPS config file - human-readable text file (tab-separated values [TSV] format) with three required and 5 optional columns:
1) *required*: LMA library name
2) *required*: barcode sequence for multiplexing - if no barcode, provide an empty string
3) *required*: complete file path to sequencing file (FASTQ or BAM)
4) *optional*: primers to exclude - comma separated list of primer names found in the second column of the primer file (**must be identical**)
5) *optional*: batch number - user-assigned batch id that groups libraries with respective input library for normalization
6) *optional*: normalization factors - one or more normalization factors found in their own additional columns (i.e., if there are three normalization factors, then there will be 8 columns present in the config file). LAMPS will multiply factors together.

&nbsp;&nbsp;&nbsp;&nbsp;*Columns 5-8 are used when normalizing 2C-ChIP libraries and not applied to 5C data

&nbsp;&nbsp;&nbsp;&nbsp;*Example and template config files for 2C-ChIP and 5C data can be found in the config_example/ folder as well as in the example  datasets

Primer file - human-readable text file (TSV format) with eight required columns:
1) *required*: genome locus name
2) *required*: primer name - for excluding primer pairs, these column values must be identical to column four's in the config file
3) *required*: primer strand - forward 'F' or reverse 'R'
4) *required*: primer sequence (**without barcode**)
5) *required*: genome assembly
6) *required*: chromosome
7) *required*: primer starting basepair
8) *required*: primer ending basepair

&nbsp;&nbsp;&nbsp;&nbsp;*Example primer files for 2C-ChIP and 5C data can be found in the example datasets

## Output

results/ - directory containing the following final output of LAMPS:
* results/interaction_matrices/*.raw.matrix - contains the frequency count of target sequences found within the sequenced library. Note - rows/cols are FWD then RVS primers resulting in the following four quadrants: F-F,F-R,R-F,R-R
* results/interaction_matrices/*log_raw.heatmap.png - heat map of the log(Frequency values) found in *.raw.matrix

&nbsp;&nbsp;&nbsp;&nbsp;*2C-ChIP specific output in the results/*:
* results/bedGraphs/raw/*.raw.bedGraph - raw frequency count of on-diagonal primer pairs in bedGraph format (easily viewable on most genome browsers)
* results/bedGraphs/norm/*.norm.bedGraph - normalized frequency count of on-diagonal primer pairs (will be input normalized if input library was found - please see LAMPS stdout to be sure, a warning will be thrown if not found)
* line_plots/< RPM or norm >/*.line_plot.png - line plots of either raw or normalized frequency counts of on-diagonal primer pairs
* results/raw_totals_vs_norm_factor.scatter.png - scatter plot of normalization factor vs. raw total counts of libraries (expected to show a linear trend - Spearman correlation provided)
* results/raw_totals_vs_norm_factor.tsv - values of normalization factor vs. raw total counts of libraries

### Secondary output

The following directories contain files specific to primer or mapping Quality Control (QC).

primer_files/ - directory containing the following files relevant to the primer pair analysis:
* *.fasta - ligated primer pair and short sequences FASTA files used as input to read aligner
* \*.matrix and *.png - primer pair similarities (BLAST bit or Bowtie 2 alignment scores between sequences) matrix and heat map. Useful for identifying problematic primers.

mapping_files/ - directory containing the following mapping and summary files:
* blastn_files/ - directory containing custom BLAST database files
* bowtie2_files/ - directory containing Bowtie 2 indices
* \*.fasta (*.fastq if conversion from BAM required) - FASTA formatted files of sequencing data used as input to read aligner
* *.mapped.tsv (BLAST) or *.mapped.sam (Bowtie2) - output of mapping read data against the ligated primer pair library
* *.mapped.log - stdout from mapping
* *.mapping_summary.bar_plot.png - plot of total and mapped reads to the ligated primer pair library
* mapping_summary.tsv - counts of total and mapped reads for libraries
* *.read_count.bar_plot.png -  bar plot breakdown of sequenced reads
* read_count_frequency_table.tsv - counts of sequenced reads breakdown
* short_read_analysis/ - contains the following files relevant to the secondary mapping of 'unmappables' to the short sequences library:
     * *.mapped.tsv (BLAST) or *.mapped.sam (Bowtie2) - output of mapping read data against the short sequences library
     * *.mapped.log - stdout from mapping
     * *_summary.bar_plot.png - plot of total and mapped reads to the short read library
     * *.short_reads_summary.tsv - counts of total and mapped reads for library

## Command line details:
Processing sequenced 2C-ChIP or 5C data:
```
LAMPS.py [-h] [--num_cpus NUM_CPUS] config primers type output

positional arguments:
  config			file path to LAMPS config file
  primers			file path to TSV file containing primer information
  type				source of LMA reads:[2C-ChIP,5C]
  output			file path to output folder	

optional arguments:
  -h, --help		show this help message and exit
  --num_cpus		set the number of cpus - default = total_num_cpus-2
  --word_size     	set the minimum required sequence length for processing (BLAST) - default = automated calculation
  --min_score		set Bowtie 2 min-score for end-to-end alignments (default = L,-0.2,-0.2)
  --no_index_build	don't re-build Bowtie 2 indices if present
```
Note - if the following error is encountered when running LAMPS with the Bowtie 2 read aligner and "-\-no_index_build" argument, the Bowtie 2 indices most likely need to be rebuilt:
```(ERR) bowtie2-align died with signal 11 (SEGV) (core dumped)```

## Testing
All software was tested on Linux Ubuntu 12.04.5 LTS (GNU/Linux 3.2.0-86-generic x86_64).

## License:
LAMPS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation.

LAMPS is distributed in the hopes that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
