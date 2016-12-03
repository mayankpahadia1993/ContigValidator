This is a validation pipeline for checking the accuracy of an assembler. 

It uses simple exact matching using a suffix tree and bwa mem and outputs the alignment percentage with a reference genome. 

It also checks how many kmers were found in the assembly as compared to the reference genome. 

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

##USAGE 

cd validationpipeline;


bash run.sh -r referenceFile.fa -s suffixTreeOutput.p -i InputFile.fa -a AlignmentResults

##OPTIONS

-r -s -i|-f -a are compulsory options. -i can be multiple


-r | --reference filename -> for passing the reference file


-s | --suffixtree filename  -> filename will store the suffixtree


-i | --input filename -> filename contains the input files to be checked


-a | --alignment filename -> filename is the file in which alignment results will be shown


-f | --file filename -> filename is a file which contains the input filenames

-suffixskip 1 -> to skip the Suffix Tree Creation. The default value is 0 and will create the Suffix Tree


-bwaskip 1 -> to skip the bwa step. The default value is 0 and it will do the bwa step


-kmer-size Integer -> The integer given here would be chosen as the kmer size. Default is 30. 


-abundance-min Integer -> The integer given here would be chosen as the abundance min. Default is 3.


##OUTPUT

The output is in the file specified with the option "-a|--alignment". The file will contain the file name and the percentage of its contigs that aligned to the reference file.

Example - 
BubblePoppedFastaFile.fa	100.0


