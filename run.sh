#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Use -h | --help for help"
    exit

fi
interactive=
inputFiles=""
inputFileArray=()
referenceGenomeArray=()
multipleGenomeFile="multipleGenomeFile.txt"
alignment="alignmentResults.txt"
suffixskip=0
cppsuffix=1
suffixsave=1
bwaskip=0
kmerskip=0
flagReference=0
flagSuffix=0
flagInput=0
flagKmer=0
clean=0
RED='\033[0;31m'
NC='\033[0m' #No Color
kmersize=30
abundancemin=1
referenceAbundanceMin=1

optionInfo="-r -s -i|-f are compulsory options. -i can be multiple \n\
-r | --reference filename -> for passing the reference file\n\
-s | --suffixtree filename  -> filename will store the suffixtree\n\
-i | --input filename -> filename contains the input files to be checked\n\
-a | --alignment filename -> filename is the file in which alignment results will be shown. Default Value is alignmentResults.txt\n\
-f | --file filename -> filename is a file which contains the input filenames\n\
-suffixskip 1 -> to skip the Suffix Tree Creation. The default value is 0 and will create the Suffix Tree\n\
-suffixsave 0 -> to not save the Suffix Tree. The default value is 1 and will save the Suffix Tree\n\
-bwaskip 1 -> to skip the bwa step. The default value is 0 and it will do the bwa step\n\
-kmerskip 1 -> to skip the kmer step. The default value is 0 and it will do the kmer step\n\
-kmer-size Integer -> The integer given here would be chosen as the kmer size. Default is 30. \n\
-abundance-min Integer -> The integer given here would be chosen as the abundance min. Default is 3. \n\
-h | --help -> for help\n\
"


while [ "$1" != "" ]; do
    case $1 in
        -r | --reference )           shift
                                referenceGenome=$1
                                referenceGenomeArray+=($1)
                                flagReference=1
                                ;;

        -s | --suffixtree )    shift
								suffixTreeOutput=$1
								flagSuffix=1
                                ;;
        -i | --input )			shift
								inputFiles="$inputFiles $1"
								inputFileArray+=($1)
								flagInput=1
								;;
		-a | --alignment )		shift
								alignment=$1
								;;
        -f | --file )			shift
								files=$1
								flagInput=1
								;;
        -suffixskip )			shift
								suffixskip=$1
								;;
		-cppsuffix )			shift
								cppsuffix=$1
								;;
		-suffixsave )			shift
								suffixsave=$1
								flagSuffix=1
								;;
		-bwaskip )				shift
								bwaskip=$1
								;;
		-kmer-size )			shift
								kmersize=$1
								flagKmer=1
								;;
		-abundance-min )		shift
								abundancemin=$1
								;;
		-kmerskip )				shift
								kmerskip=$1
								flagKmer=1
								;;
		-clean )				shift
								clean=$1
								;;
        -h | --help )           printf -- "$optionInfo"
                                exit
                                ;;
        * )                     printf -- "$optionInfo"
                                exit 
    esac
    shift
done

# echo "here"
if [ "$flagReference" = 0 ]; then
	echo "You didn't set the reference flag (-r | --reference)"
	printf -- "$optionInfo"
	exit
fi

if [ "$flagSuffix" = 0 ]; then
	echo "You didn't set the suffix tree flag (-s | --suffixtree)"
	printf -- "$optionInfo"
	exit
fi

if [ "$flagInput" = 0 ]; then
	echo "You didn't set the input flag (-i | --input) or the file flag(-f | --file). You need to set one of them."
	printf -- "$optionInfo"
	exit
fi

if [ "$flagKmer" = 0 ]; then
	echo "You didn't set the kmerskip flag (-kmerskip) or the kmer-size flag(-kmer-size). You need to set one of them."
	printf -- "$optionInfo"
	exit
fi

if [ "$abundancemin" -le 0 ]; then
	echo "Abundance min can't be 0 or negative"
	exit
fi

if [ "$kmersize" -le 0 ]; then
	echo "Kmer Size can't be 0 or negative"
	exit
fi

if [ -n "$files" ]; then
    # echo "You have set the -f | --file option. Hence we will check the files mentioned in $files"
    inputFiles=""
    inputFileArray=()
    while read line; do    
    inputFiles+="$inputFiles $line"
    inputFileArray+=($line)
	done < $files
fi

## For multiple genomes, we concatenate all of them to one file. 

## to empty the contents of the file
> $multipleGenomeFile

##concatenate every file to multipleGenomeFile
for i in "${referenceGenomeArray[@]}"
do
	cat $i >> $multipleGenomeFile
done

##the first reference file given is assumed as the main reference file without any mutations or substitutions
referenceGenome=${referenceGenomeArray[0]}

## Suffix Tree for simple exact matching
if [ "$cppsuffix" = 1 ]; then
	echo "Running C++ Suffix Tree"
	g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib src/program.cpp -o src/program -lsdsl -ldivsufsort -ldivsufsort64
	## Adding filenames to the $alignment output file
	echo "Filename" > $alignment
	for i in "${inputFileArray[@]}"
	do
		echo $i >> "$alignment"
	done


	## Check of input files with simple exact matching

	echo "Checking the input files against the tree"
	# python src/finalAlign.py $referenceGenome tempSuffix.txt $inputFiles > tempout.txt 2> tempError.txt
	# src/program $referenceGenome tempSuffix.txt $inputFiles > tempout.txt 2> tempError.txt
	# src/program $referenceGenome tempSuffix.txt $inputFiles > tempout.txt 2> tempError.txt

	##Performing it with multipleGenomeFile
	src/program $multipleGenomeFile tempSuffix.txt $inputFiles > tempout.txt 2> tempError.txt
	
	# time ~/medvedevGroup/playWithSDSL/program ../ecoliReference.fa ecoliReadsNew150-0-1BubblePoppedFastaFile.fa log.txt

	if [ "$?" = 0 ]; then
		cat tempout.txt
		paste $alignment tempSuffix.txt > tempout.txt
		cat tempout.txt > $alignment
		rm -f tempSuffix.txt
	else
		printf "${RED}ERROR - "
		# echo -e "I ${RED}love${NC}"
		cat tempError.txt
		printf "${NC}"
	fi
	echo "Exact Alignment done"
else
	if [ "$suffixsave" = 1 ]; then
		if [ "$suffixskip" = 0 ]; then
			echo "Creating the Suffix Tree"
			# python src/createSuffixTree.py $referenceGenome $suffixTreeOutput > tempout.txt 2> tempError.txt

			##Performing it with multipleGenomeFile
			python src/createSuffixTree.py $multipleGenomeFile $suffixTreeOutput > tempout.txt 2> tempError.txt

			if [ "$?" = 0 ]; then
				cat tempout.txt
			else
				printf "${RED}ERROR - "
				# echo -e "I ${RED}love${NC}"
				cat tempError.txt
				printf "${NC}"
			fi
		else
			echo "Skipping Suffix Tree creation and using the suffix tree in $suffixTreeOutput"
		fi


		## Adding filenames to the $alignment output file
		echo "FILENAME" > $alignment
		for i in "${inputFileArray[@]}"
		do
			echo $i >> "$alignment"
		done


		## Check of input files with simple exact matching

		echo "Checking the input files against the tree"
		python src/checkWithSuffixTree.py $suffixTreeOutput tempSuffix.txt $inputFiles > tempout.txt 2> tempError.txt

		if [ "$?" = 0 ]; then
			cat tempout.txt
			paste $alignment tempSuffix.txt > tempout.txt
			cat tempout.txt > $alignment
			rm -f tempSuffix.txt
		else
			printf "${RED}ERROR - "
			# echo -e "I ${RED}love${NC}"
			cat tempError.txt
			printf "${NC}"
		fi

	else
		echo "Not saving Suffix Tree"
		## Adding filenames to the $alignment output file
		echo "Filename" > $alignment
		for i in "${inputFileArray[@]}"
		do
			echo $i >> "$alignment"
		done


		## Check of input files with simple exact matching

		echo "Checking the input files against the tree"
		# python src/finalAlign.py $referenceGenome tempSuffix.txt $inputFiles > tempout.txt 2> tempError.txt

		##Performing it with multipleGenomeFile
		python src/finalAlign.py $multipleGenomeFile tempSuffix.txt $inputFiles > tempout.txt 2> tempError.txt

		if [ "$?" = 0 ]; then
			cat tempout.txt
			paste $alignment tempSuffix.txt > tempout.txt
			cat tempout.txt > $alignment
			rm -f tempSuffix.txt
		else
			printf "${RED}ERROR - "
			# echo -e "I ${RED}love${NC}"
			cat tempError.txt
			printf "${NC}"
		fi
	fi
fi
	





## BWA

if [ "$bwaskip" = 0 ]; then
	echo "BWA - "
	# Index the reference file

	bwa index $referenceGenome > tempout.txt 2> tempError.txt

	if [ "$?" = 0 ]; then
		cat tempout.txt
	else
		printf "${RED}ERROR - "
		# echo -e "I ${RED}love${NC}"
		cat tempError.txt
		printf "${NC}"
	fi

	# align the input files with the reference file
	bwaOutput=""
	echo "%align" > tempBwaOutput.txt
	for i in "${inputFileArray[@]}"
	do
		output="$i.bwa.bam"
		bwa mem $referenceGenome $i | samtools sort > $output 2> tempError.txt
		if [ "$?" -gt 0 ]; then
			printf "${RED}ERROR - "
			# echo -e "I ${RED}love${NC}"
			cat tempError.txt
			printf "${NC}"
		fi
		samtools index $output 2> tempError.txt
		if [ "$?" -gt 0 ]; then
			printf "${RED}ERROR - "
			# echo -e "I ${RED}love${NC}"
			cat tempError.txt
			printf "${NC}"
		fi
		samtools flagstat $output | grep mapped | cut -f 5 -d " "| cut -f 2 -d "(" | head -1 >> tempBwaOutput.txt 2> tempError.txt
		if [ "$?" -gt 0 ]; then
			printf "${RED}ERROR - "
			# echo -e "I ${RED}love${NC}"
			cat tempError.txt
			printf "${NC}"
		fi
		# echo $tempBwaOutput >> tempBwaOutput.txt
		# bwaOutput+="$i"
		# bwaOutput+="\t"
		# bwaOutput+=$tempBwaOutput
		# bwaOutput+="\n"

	done

	# echo -e -n "$bwaOutput"

	paste $alignment tempBwaOutput.txt > tempout.txt
	cat tempout.txt > $alignment
	rm -f tempBwaOutput.txt

else
	echo "Skipping BWA"
fi

if [ "$kmerskip" = 0 ]; then

	echo "Kmer size = $kmersize"
	echo "Abundance Min = $abundancemin"

	## dsk on reference 
	h5file="$referenceGenome.h5"
	# dsk -file $referenceGenome -kmer-size $kmersize -abundance-min $referenceAbundanceMin -out $h5file 2> tempError.txt

	##Performing it with multipleGenomeFile
	dsk -file $multipleGenomeFile -kmer-size $kmersize -abundance-min $referenceAbundanceMin -out $h5file 2> tempError.txt


	if [ "$?" -gt 0 ]; then
			printf "${RED}ERROR - "
			# echo -e "I ${RED}love${NC}"
			cat tempError.txt
			printf "${NC}"
	fi
	echo "dsk done"

	kmercountfile="$referenceGenome.kmercount"
	dsk2ascii -file $h5file -out $kmercountfile 2> tempError.txt
	if [ "$?" -gt 0 ]; then
			printf "${RED}ERROR - "
			# echo -e "I ${RED}love${NC}"
			cat tempError.txt
			printf "${NC}"
	fi
	echo "dsk2ascii done"
	if [ "$clean" = 0 ]; then
		rm -f $h5file
	fi

	## dsk on files

	for i in "${inputFileArray[@]}"
	do
		h5file="$i.h5"
		dsk -file $i -kmer-size $kmersize -abundance-min $abundancemin -out $h5file 2> tempError.txt
		if [ "$?" -gt 0 ]; then
				printf "${RED}ERROR - "
				# echo -e "I ${RED}love${NC}"
				cat tempError.txt
				printf "${NC}"
		fi
		h5file="$i.h5"
		kmercountfile="$i.kmercount"
		dsk2ascii -file $h5file -out $kmercountfile 2> tempError.txt
		if [ "$?" -gt 0 ]; then
				printf "${RED}ERROR - "
				# echo -e "I ${RED}love${NC}"
				cat tempError.txt
				printf "${NC}"
		fi
		if [ "$clean" = 0 ]; then
			rm -f $h5file
		fi
	done

	##Setting up the input for findCommonKmers.py file

	commonKmerInputFile=""
	inputkmercountfile=""
	for i in "${inputFileArray[@]}"
	do
		file1="$i.commonKmers12"
		file2="$i.commonKmers21"
		inputkmercountfile+=" $i.kmercount"
		commonKmerInputFile+=" $i.kmercount $file1 $file2"
	done

	## Calling findCommonKmers.py file

	tempKmerOut="tempKmerOut.txt"

	# python src/findCommonKmers.py $tempKmerOut "$referenceGenome.kmercount" $commonKmerInputFile > tempout.txt 2> tempError.txt
	## Change to C++ code for findCommonKmers.cpp here

	
	g++ -std=c++0x -o src/findCommonKmers src/findCommonKmers.cpp
	src/findCommonKmers "$referenceGenome.kmercount" $tempKmerOut $inputkmercountfile > tempout.txt 2> tempError.txt


	if [ "$?" = 0 ]; then
		cat tempout.txt
		paste $alignment $tempKmerOut > tempout.txt
		cat tempout.txt > $alignment
		rm -f $tempKmerOut
	else
		printf "${RED}ERROR - "
		# echo -e "I ${RED}love${NC}"
		cat tempError.txt
		printf "${NC}"
	fi

	rm -f tempout.txt
	rm -f tempError.txt

else
	echo "Skipping Kmer step"
fi

if [ "$clean" = 0 ]; then
	echo "In clean"
	rm -f $kmercountfile
	rm -f $inputkmercountfile
	rm -f $referenceGenome.*
	# rm -f $multipleGenomeFile
fi


# python src/findCommonKmers.py $tempKmerOut "$referenceGenome.kmercount" $commonKmerInputFile
