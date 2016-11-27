This is a validation pipeline for checking the accuracy of an assembler. 

It uses simple exact matching using a suffix tree and bwa mem and outputs the alignment percentage with a reference genome. 

It also checks how many kmers were found in the assembly as compared to the reference genome. 

To get the kmer counts, we use dsk. 

Suffix tree is from - https://github.com/kvh/Python-Suffix-Tree

We also use bwa mem to align and get the percentage. 
