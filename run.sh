#!/bin/bash

tmpoutfile=tmp_out.txt
tmperrorfile=tmp_err.txt
tmpsuffixfile=temp_suffix.txt
tmpbwafile=tempBwaOutput.txt

SCRIPT_PATH="$(cd "$(dirname $0)"; pwd)"

function sane_quit() {
    rm -f $tmpoutfile $tmperrorfile $mutipleGenomFile $tmpsuffixfile $tmpbwafile
    exit
}

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
cat ${referenceGenomeArray} > $multipleGenomeFile

##the first reference file given is assumed as the main reference file without any mutations or substitutions
referenceGenome=${referenceGenomeArray[0]}

## Suffix Tree for simple exact matching
echo "Running C++ Suffix Tree"
cd src
make 
if [ $? -ne 0 ]; then
    echo "Compilation failed.  Please check that sdsl is in your include path and in your library path."
    cd ..
    sane_quit
fi
cd ..

## Adding filenames to the $alignment output file
echo "Filename" > $alignment

for i in "${inputFileArray[@]}"
do
    echo $i >> "$alignment"
done

## Check of input files with simple exact matching
echo "Checking the input files against the tree"

##Performing it with multipleGenomeFile
src/program $multipleGenomeFile $tmpsuffixfile $inputFiles > $tmpoutfile 2> $tmperrorfile

if [ $? = 0 ]; then
    cat $tmpoutfile
    paste $alignment $tmpsuffixfile > $tmpoutfile
    cat $tmpoutfile > $alignment
    rm -f $tmpsuffixfile
    rm -f $tmpoutfile
    rm -f $tmperrorfile
else
    printf "${RED}ERROR - "
    cat $tmperrorfile
    printf "${NC}"
    sane_quit
fi
echo "Exact Alignment done"

## BWA
if [ "$bwaskip" = 0 ]; then
    echo "BWA - "
    # Index the reference file
    bwa index $referenceGenome > $tmpoutfile 2> $tmperrorfile
    
    if [ "$?" = 0 ]; then
	cat $tmpoutfile
    else
	printf "${RED}ERROR - "
	cat $tmperrorfile
	printf "${NC}"
        sane_quit
    fi
    
    # align the input files with the reference file
    bwaOutput=""
    echo "%align" > $tmpbwafile
    for i in "${inputFileArray[@]}"
    do
	output="$i.bwa.bam"
	bwa mem -t 16 $referenceGenome $i | samtools sort > $output 2> $tmperrorfile
	if [ "$?" -gt 0 ]; then
	    printf "${RED}ERROR - "
	    cat $tmperrorfile
	    printf "${NC}"
            sane_quit
	fi
	samtools index $output 2> $tmperrorfile
	if [ "$?" -gt 0 ]; then
	    printf "${RED}ERROR - "
	    cat $tmperrorfile
	    printf "${NC}"
            sane_quit
	fi
	samtools flagstat $output | grep mapped | cut -f 5 -d " "| cut -f 2 -d "(" | head -1 >> $tmpbwafile 2> $tmperrorfile
	if [ "$?" -gt 0 ]; then
	    printf "${RED}ERROR - "
	    # echo -e "I ${RED}love${NC}"
	    cat $tmperrorfile
	    printf "${NC}"
	fi
    done
    
    paste $alignment $tmpbwafile > $tmpoutfile
    cat $tmpoutfile > $alignment
    rm -f $tmpbwafile
else
    echo "Skipping BWA"
fi

# USE KMC FOR BETTER PERFOMANCE
cp $alignment old-alignment.txt

cd ${SCRIPT_PATH}/src
make
cd ${SCRIPT_PATH}

tmpKmerOut="/tmp/tmpkmerout.txt"
if [ "$kmerskip" = 0 ]; then
    echo "Kmer size = $kmersize"
    echo "Abundance Min = $abundancemin"

    ${SCRIPT_PATH}/src/KMC/bin/kmc -t4 -ci1 -k$kmersize -fm $multipleGenomeFile mg /tmp

    echo -e "recall\tprecision" > $tmpKmerOut
    for i in "${inputFileArray[@]}"
    do
        ${SCRIPT_PATH}/src/KMC/bin/kmc -t4 -ci$abundancemin -k$kmersize -fm $i $i.kmc /tmp
        ${SCRIPT_PATH}/src/KMC/bin/kmc_tools simple mg -ci1 $i.kmc -ci1 intersect $i.tp -ci1
        ${SCRIPT_PATH}/src/KMC/bin/kmc_tools simple mg -ci1 $i.kmc -ci1 kmers_subtract $i.fn -ci1
        ${SCRIPT_PATH}/src/KMC/bin/kmc_tools simple mg -ci1 $i.kmc -ci1 reverse_kmers_subtract $i.fp -ci1
        tp=$(./src/count_kmers_kmc $i.tp)
        fn=$(./src/count_kmers_kmc $i.fn)
        fp=$(./src/count_kmers_kmc $i.fp)
        precision=$(bc -l <<< "scale=4; $tp * 100 / ($tp + $fp)")
        recall=$(bc -l <<< "scale=4; $tp * 100 / ($tp + $fn)")
        echo -e "$recall%\t$precision%" >> $tmpKmerOut
    done
    paste $alignment $tmpKmerOut > $tmpoutfile
    cat $tmpoutfile > $alignment
fi

if [ "$clean" = 0 ]; then
    echo "In clean"
    rm -f $referenceGenome.*
    rm -f *.kmc_pre *.kmc_suf
fi

sane_quit
