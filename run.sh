#!/bin/bash

if [ $# -lt 4 ]; then
    echo "Your command line contains less arguments"
    echo "-r -s -i|-f -a are compulsory options. -i can be multiple"
    echo "Use -r | --reference filename -> for passing the reference file"
    echo "Use -s | --suffixtree filename  -> filename will store the suffixtree"
    echo "Use -i | --input filename -> filename contains the input files to be checked"
    echo "Use -a | --alignment filename -> filename is the file in which alignment results will be shown"
    echo "Use -f | --file filename -> filename is a file which contains the input filenames"

fi
interactive=
inputFiles=()


while [ "$1" != "" ]; do
    case $1 in
        -r | --reference )           shift
                                referenceGenome=$1
                                ;;

        -s | --suffixtree )    shift
								suffixTreeOutput=$1
                                ;;
        -i | --input )			shift
								inputFiles+=($1)
								;;
		-a | --alignment )		shift
								alignment=$1
								;;
        -f | --file )			shift
								files=$1
								;;
        -h | --help )           echo "-r -s -i|-f -a are compulsory options. -i can be multiple"
							    echo "Use -r | --reference filename -> for passing the reference file"
							    echo "Use -s | --suffixtree filename  -> filename will store the suffixtree"
							    echo "Use -i | --input filename -> filename contains the input files to be checked"
							    echo "Use -a | --alignment filename -> filename is the file in which alignment results will be shown"
							    echo "Use -f | --file filename -> filename is a file which contains the input filenames"
                                exit
                                ;;
        * )                     echo "-r -s -i|-f -a are compulsory options. -i can be multiple"
							    echo "Use -r | --reference filename -> for passing the reference file"
							    echo "Use -s | --suffixtree filename  -> filename will store the suffixtree"
							    echo "Use -i | --input filename -> filename contains the input files to be checked"
							    echo "Use -a | --alignment filename -> filename is the file in which alignment results will be shown"
							    echo "Use -f | --file filename -> filename is a file which contains the input filenames"
                                exit 
    esac
    shift
done

echo "here"



echo "Creating the Suffix Tree"

python src/createSuffixTree.py $referenceGenome $suffixTreeOutput

if [ -n "$files" ]; then
    echo "You have set the -f | --file option. Hence we will check the files mentioned in $files"
    inputFiles=()
    while read line; do    
    inputFiles+=($line)
	done < $files
fi


o="Output"
for i in "${inputFiles[@]}"
do
	echo $i
	output=$i$o
	echo "Checking the input files against the tree"
	echo $output
	python src/checkWithSuffixTree.py $suffixTreeOutput $i $output $alignment
done






