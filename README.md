#ContigValidator

This is a pipeline for validating the accuracy of a set of contigs with respect to a known referene genome. 
The set of contigs can be just the read themselves, or unitigs, or any other type of contig.
Given a fasta reference genome and a file containing a set of contigs, ContigValidator reports:

1. The percentage of contigs that occur exactly in the reference
2. The percentage of contigs that align to the reference.
3. Identifies k-mers shared between the contigs and the reference and reports the corresponding percentages.


<!-- This is a validation pipeline for checking the accuracy of an assembler. 

The inputs to the pipeline are fasta files. A fasta file for the reference genome and one or more fasta files for contigs/reads.

This pipeline does the following things - 

1. Uses simple exact matching to get the alignment percent of a read/contig with the reference
2. Uses BWA to find the alignment of read/contig with the reference
3. Compares the percentage of kmers common between the reads/contig and reference
--->


#Release Date
###TBD

#Prerequisites

1. Python 2.7 - Should be present in the PATH variable


2. BWA - version used bwa-0.7.15-r1140, Should be present in the PATH variable

		# Install using homebrew
		brew install bwa

3. DSK - version used 2.1.0, Should be present in the PATH variable

		# Install using homebrew
		brew install dsk

4. SAMTOOLS - version used 1.3.1, Should be present in the PATH variable

		# Install using homebrew
		brew install samtools

#Installation

To download the program, you can clone the repository from github, either using SSH: 

	
	git clone git@github.com:mayankpahadia1993/ContigValidator.git
	
or HTTPS
	
	git clone https://github.com/mayankpahadia1993/ContigValidator.git

No installation is needed and the program is ready to run.

#Testing

To make sure that the program was installed correctly and works as intended, you can test it:

	cd validationpipeline/test/;
	bash test.sh;

The output should  report the message "TEST PASSED SUCCESSFULLY" at the end to indicate that the installation is working properly. Otherwise, an error message should appear indicating the source of the error. 


#Example use case

Lets say you would like to test the accuracy of your new contig assembler. 
You took the Ebola reference sequence, stored as reference.fa, and 
simulated a dataset of reads, stored as reads.fa. 
Initially, you can validate the accuracy of the reads using the command

	bash run.sh -r reference.fa -s ebola_suffix_tree -i reads.fa -a results.reads.txt -kmer-size 30 

After running this command, we can check the results ```cat results.reads.txt``` and get:
	
	Filename	%exact	%align	recall	precision
	reads.fa	22.136	100.00%	100.0%	1.27982540829%

Next, we assemble the reads using your new contig assembler. 
You want to try three different error correction strategies, so you get three different contig files: 
`contig1.fa`, `contig2.fa`, `contig3.fa`.
You can run the validation pipeline as follows:
	
	bash run.sh -r reference.fa -s ebola_suffix_tree -suffixskip 1 -a results.assembly.txt -kmer-size 30 -i contig1.fa -i contig2.fa -i contig3.fa

After running this command, we can check the results ```cat results.assembly.txt``` and get:
	
	Filename	%exact	%align	recall	precision
	contig1.fa	
	contig2.fa	
	contig3.fa	
	


#Detailed options


`bash run.sh  -r <filename> [-s <filename>] [-suffixsave <0|1>] [-i <filename>]  [-f <filename>]  [-a <filename>] [-suffixskip <0|1>] [-bwaskip <0|1>] [-kmerskip <0|1>] [-kmer-size <int>] [-abundance-min <int>]`


Where: 

`-r <filename>`, `--reference <filename>`: the reference genome. This should be a fasta file.


`-s <filename>`, `--suffixtree <filename>`: The location of the suffix tree. If the -suffixskip=0 option is used, this is the location of the output file containing the suffix tree. If the -suffixskip=1 option is used, this is the location of the input file which contains the suffix tree. If -suffixsave=0 option is used then this option neeed not be provided.


`-suffixsave <0|1>`: ContigValidator needs to build a suffix tree for the reference file before doing alignment. By default (value 1), the suffix tree is created on the fly and written to the file specified by the -s option so that it can be reused in future runs. If, on the other hand, you don't need to save the suffixtree into a file, you can use the -suffixsave=0 option to not save it. 


`-i <filename>`, `--input <filename>`: The location of the contig file. This should be a fasta file, with each entry corresponding to a contig. This option must be specified if the -f option is not specified. The -i option may be used multiple times to specify multiple input files.


`-f <filename>`, `--file <filename>`: A file which contains the names of contig files, one per line. This option can be used in place of the -i option when there are multiple files. This option must be specified if the -i option is not specified.


`-a <filename>`, `--alignment <filename>`: The location of the output file with the results of the analysis (default: alignmentResults.txt)


`-suffixskip <0|1>`: ContigValidator needs to build a suffix tree for the reference file before doing alignment. By default (value 0), the suffix tree is created on the fly and written to the file specified by the -s option so that it can be reused in future runs. If, on the other hand, a suffix tree file is already available from a previous run, then you can set suffixskip to 1 in order to skip the creation of the suffix tree and instead load it from the fly specified by the -s option.


`-bwaskip <0|1>`: ContigValidator uses BWA to align the contigs to the reference. To skip this step, set the value of this option to 1.  (Default: 0)



`-kmerskip <0|1>`: ContigValidator uses kmer counts to figure out recall and precision. To skip this step, set the value of this option to 1. This option must be specified if the -kmer-size option has not been specified.


`-kmer-size <int>` The size of the kmers used. This option must be specified if the -kmerskip option has not been specified. 


`-abundance-min <int>`: Kmers that occur in the contigs less than abundance-min times are treated as not present in the contig file.  (Default: 3) 


### Help Information:

use `-h/--help` for detailed help message.


#Interpreting the results

These are the output files that are generated by the pipeline - 

1. AlignmentResults
2. *.commonKmers12
3. *.commonKmers21
4. SuffixTreeOutput file
5. BAM File
6. *.exact


###1. AlignmentResults 

The output is in the file specified with the option "-a". Each row represents an input file, as given by -i or -f in the command line. For each row, there are up to 5 tab-delimited columns:


- Filename: the name of the input contig file.


- %exact: the percentage of contigs that are exactly present (as substrings) in either the reference or the reverse complement of the reference.


- %align: the percentage of contigs for which there exists an alignment to the reference, using BWA MEM. This percentage is taken from the "mapped percentage" printed in the bwa output.


- %recall: the percentage of k-mers from the reference file that are also abundandtly present (exactly) in the contig input file. A kmer is abundandatly present if it occurs at least the number of times specified by -abundance-min. 


- %precision: the percentage of k-mers abundandtly present in the contig input file that are also present (exactly) in the reference file.


Example - 


FILENAME	%exact	%align	recall	precision

test/reads.fa	22.136	100.00%	100.0%	1.27982540829%


###2. *.commonKmers12 
These files are generated for all the contig input files. In this file, there is one row per abundant kmer from the contig input file. Every row contains two tab-delimited columns. The first column is the kmer and the second column is a binary number. '1' represents the kmer is present in the reference file and '0' represents it is not.


###3. *.commonKmers21
These files are generated for all the contig input files. In this file, there is one row per kmer from the reference file. Every row contains two tab-delimited columns. The first column is the kmer and the second column is a binary number. '1' represents the kmer is abundantly present in the contig input file and '0' represents it is not. 


###5. *.bwa.bam
These files are generated for all the contig inputfile, if the bwaskip option is not set. They store the output of bwa-mem.

###6. *.exact

This file is generated for all the input contig files. In the file, there is one row per contig. Each row has two tab-delimited columns. The first column is the contig id from the contig fasta file. The second column is a binary number. '1' represents that the contig is present exactly in the reference or its reverse complement. '0' represents that it is not. 


#Authors
- Mayank Pahadia (The Pennsylvania State University)
- Paul Medvedev (The Pennsylvania State University)


# Contact

- You  can report bugs or suggest features using issue tracker at GitHub [https://github.com/mayankpahadia1993/validationpipeline](https://github.com/mayankpahadia1993/validationpipeline)

- You can contact Mayank Pahadia directly at mqp5492@cse.psu.edu


