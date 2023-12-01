# Hossain_etal_2024
Custom Python scripts used in Hossain et al., 2024 for analysis of next-generation sequencing (NGS) data and NCBI protein BLAST XML alignment data:

NGS_phage_DNA_infected_cells.py takes as input a FASTA file with a phage reference genome sequence as well as two NGS Illumina FASTQ files with total DNA reads from bacterial cells infected with phage, under control or immune targeting conditions, to map the base pair coverage across the entire phage genome normalized to the total DNA reads for a particular sample. The output is a graph of the data as well as two CSV files for each FASTQ input. One CSV file records all the DNA reads that mapped to the phage genome with their start and end nucleotide position numbers indicated. The other CSV file records the normalized read count for each base per position within the phage genome.

Also provided is a test client for the Python program above, titled pWEBTNC_pBrig1_T4_test_client.py, which runs the program above on FASTQ reads for pWEB-TNC and pBrig1, control and targeting conditions, respectively, aligned to a provided reference genome for wild-type T4 phage (T4_WT.fa). FASTQ reads represent total DNA extracted from E. coli cells infected with wild-type T4 phage for 8 minutes. For pWEB-TNC (sample S1), reads are from cells carrying the control empty cosmid, pWEB-TNC. For pBrig1 (sample S2), reads are from cells expressing the Brig1 DNA glycosylase from the cosmid pBrig1. Two reference phage genome FASTA files are provided in this repository, T4_WT.fa and T4_escaper1.fa. FASTQ files can be found at the NCBI Sequence Read Archive (SRA) under NCBI BioProject ID PRJNA1045662.

Gene_Neighborhood_Analysis.py takes as input a protein BLAST or PSI-BLAST XML file and returns as output a TSV file with n genes upstream and n genes downstream of each protein accession hit in the XML file by interfacing with the Entrez database and locating the position of each BLAST hit within its associated nucleotide entry. The PSI-BLAST alignment XML file (NZUWTKTR016-Alignment.xml) used for Brig1 gene neighborhood analysis is also provided in this repository.


HARDWARE REQUIREMENTS

Both scripts require only a standard computer with enough RAM to support the in-memory operations.


SOFTWARE REQUIREMENTS

The code is supported for standard operating systems with functionality for Python. Scripts have been tested on the following systems: macOS Sonoma version 14.1.2. Scripts were written and tested using the software PyCharm CE (PyCharm 2023.2.5 Community Edition). Code was written and tested using Python 3.1 (2).


PYTHON DEPENDENCIES
The Python scripts in this repository mainly depend on the following:
numpy
biopython
pandas
matplotlib


RUNTIME AND REPRODUCIBILITY

The runtime for the pWEBTNC_pBrig1_T4_test_client.py, which uses only forward NGS reads from a paired-end 2 x 75 bp Illumina MiSeq run, is ~2 hours using PyCharm CE (PyCharm 2023.2.5 Community Edition) on a 13-inch 2019 MacBook Pro with a 2.8 GHz Quad-Core Intel Core i7 processor and 16 GB 2133 MHz LPDDR3 running on macOS Sonoma version 14.1.2. The FASTA file for this test client is "T4_WT.fa" available in this repository. 


LICENSE
This project is covered under the Apache 2.0 License.

 Copyright [2023] [Amer A. Hossain, Christian F. Baca, Luciano A. Marraffini]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.








