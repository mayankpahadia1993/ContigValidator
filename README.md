This is a validation pipeline for checking the accuracy of an assembler. 

This pipeline does the following things - 

1. [Uses simple exact matching to get the alignment percent of a read/contig with the reference]
2. [Uses BWA to find the alignment of read/contig with the reference]
3. [Compares the percentage of kmers common between the reads/contig and reference]

For Simple Exact Matching, we use a suffix tree and create a suffix tree of the reference and check if the contigs/reads in the file are present in the suffix tree or not.

To get the kmer counts, we use dsk - https://github.com/GATB/dsk

Suffix tree is from - https://github.com/kvh/Python-Suffix-Tree

We also use bwa mem to align and get the percentage - https://github.com/lh3/bwa

##PREREQUISITES -

Python 2.7


BWA - version used bwa-0.7.15-r1140, Should be present in the PATH

Install using homebrew - brew install bwa

DSK - version used 2.1.0, Should be present in the PATH

Install using homebrew - brew install dsk

SAMTOOLS - version used 1.3.1, Should be present in the PATH

Install using homebrew - brew install samtools

##INSTALLATION

To use the program

SSH - 
			git clone git@github.com:mayankpahadia1993/validationpipeline.git


HTTPS - 
			git clone https://github.com/mayankpahadia1993/validationpipeline.git

##TESTING
			
			cd validationpipeline/test/;
			bash test.sh;

It should print "EVERYTHING IS WELL" at the end. Or otherwise it will give an error.


##USAGE 


			bash run.sh -r referenceFile.fa -s suffixTreeOutput.p -i InputFile.fa -a AlignmentResults

##OPTIONS

-r -s -i|-f are compulsory options. -i can be multiple


-r | --reference filename -> for passing the reference file


-s | --suffixtree filename  -> filename will store the suffixtree


-i | --input filename -> filename contains the input file to be checked


-a | --alignment filename -> filename is the file in which alignment results will be shown. Default Value is alignmentResults.txt


-f | --file filename -> filename is a file which contains the input filenames


-suffixskip 1 -> to skip the Suffix Tree Creation. The default value is 0 and will create the Suffix Tree


-bwaskip 1 -> to skip the bwa step. The default value is 0 and it will do the bwa step


-kmer-size Integer -> The integer given here would be chosen as the kmer size. Default is 30. 


-abundance-min Integer -> The integer given here would be chosen as the abundance min. Default is 3.


##INPUT

The mandatory input to the program is the reference file, suffix tree output file and the files containing reads and contigs.

The reference file and read files must be in the fasta format. The reads here usually mean contigs or unitigs which we get after running an assembler.

##OUTPUT

The output is in the file specified with the option "-a|--alignment". The file will contain a maximum of 5 columns all separated by tabs.

Example - 
FILENAME	Align Percentage	BWA Percentage	Percentage of Kmers from Reference File in Input File	Percentage of Kmers from Input File in Reference File
test/reads.fa	22.136	100.00%	100.0%	1.27982540829%

The output has 5 columns - Filename, Align Percentage, BWA Percentage, Percentage of Kmers from Reference File in Input File and Percentage of Kmers from Input File in Reference File

1. [Filename] - This column has the name of the input files(reads/contigs).
2. [Align Percentage] - This column has the align percentage calculated using simple exact matching.
3. [BWA Percentage] - This column has the align percentage calculated using BWA MEM.
4. [Percentage of Kmers from Reference File in Input File] - This column has the percentage of kmers common from the reference file in the input file.
5. [Percentage of Kmers from Input File in Reference File] - This column has the percentage of kmers common from the input file in the reference file.

