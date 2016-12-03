#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Use -h | --help for help"
    exit

fi
interactive=
inputFiles=""
skip=0


while [ "$1" != "" ]; do
    case $1 in
        -r | --reference )           shift
                                referenceGenome=$1
                                ;;

        -s | --suffixtree )    shift
								suffixTreeOutput=$1
                                ;;
        -i | --input )			shift
								inputFiles="$inputFiles $1"
								;;
		-a | --alignment )		shift
								alignment=$1
								;;
        -f | --file )			shift
								files=$1
								;;
        -skip )					shift
								skip=$1
								;;
        -h | --help )           echo "-r -s -i|-f -a are compulsory options. -i can be multiple"
							    echo "Use -r | --reference filename -> for passing the reference file"
							    echo "Use -s | --suffixtree filename  -> filename will store the suffixtree"
							    echo "Use -i | --input filename -> filename contains the input files to be checked"
							    echo "Use -a | --alignment filename -> filename is the file in which alignment results will be shown"
							    echo "Use -f | --file filename -> filename is a file which contains the input filenames"
							    echo "Use -skip 1 -> to skip the Suffix Tree Creation. The default value is 0 and will create the Suffix Tree"
							    echo "Use -h | --help -> for help"
                                exit
                                ;;
        * )                     echo "-r -s -i|-f -a are compulsory options. -i can be multiple"
							    echo "Use -r | --reference filename -> for passing the reference file"
							    echo "Use -s | --suffixtree filename  -> filename will store the suffixtree"
							    echo "Use -i | --input filename -> filename contains the input files to be checked"
							    echo "Use -a | --alignment filename -> filename is the file in which alignment results will be shown"
							    echo "Use -f | --file filename -> filename is a file which contains the input filenames"
							    echo "Use -skip 1 -> to skip the Suffix Tree Creation. The default value is 0 and will create the Suffix Tree"
							    echo "Use -h | --help -> for help"
                                exit 
    esac
    shift
done

# echo "here"

if [ -n "$files" ]; then
    echo "You have set the -f | --file option. Hence we will check the files mentioned in $files"
    inputFiles=""
    while read line; do    
    inputFiles+="$inputFiles $line"
	done < $files
fi


if [ "$skip" = 0 ]; then
	echo "Creating the Suffix Tree"
	python src/createSuffixTree.py $referenceGenome $suffixTreeOutput
else
	echo "Skipping Suffix Tree creation and using the suffix tree in $suffixTreeOutput"
fi


echo "Checking the input files against the tree"
python src/checkWithSuffixTree.py $suffixTreeOutput $alignment $inputFiles







