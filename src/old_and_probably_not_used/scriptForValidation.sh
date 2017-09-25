echo "Alignment of Reference to bcalm output"

python ~/medvedevGroup/bubblepopping/src/align.py ../ecoliReference.fa ecoliReadsNew150-0-1.unitigs.fa bcalmResultsSuffix.txt > bcalmConsoleOutputResultsAlignment.txt

echo "Alignment of Reference to bubblepopping output"

python ~/medvedevGroup/bubblepopping/src/align.py ../ecoliReference.fa ecoliReadsNew150-0-1BubblePoppedFastaFile.fa bubblePoppedResultsSuffix.txt > bubblepoppingConsoleOutputResultsAlignment.txt

echo "Kmer Counting on Bcalm"

~/medvedevGroup/dsk/build/bin/dsk -file ecoliReadsNew150-0-1.unitigs.fa -kmer-size 30 -abundance-min 1

~/medvedevGroup/dsk/build/bin/dsk2ascii -file ecoliReadsNew150-0-1.unitigs.h5 -out ecoliReadsNew150-0-1BcalmKmerCount.txt

echo "Kmer Counting on bubblepopping"

~/medvedevGroup/dsk/build/bin/dsk -file ecoliReadsNew150-0-1BubblePoppedFastaFile.fa -kmer-size 30 -abundance-min 1

~/medvedevGroup/dsk/build/bin/dsk2ascii -file ecoliReadsNew150-0-1BubblePoppedFastaFile.h5 -out ecoliReadsNew150-0-1BubblePoppedKmerCount.txt

echo "Kmer Counting on reference"

~/medvedevGroup/dsk/build/bin/dsk -file ../ecoliReference.fa -kmer-size 30 -abundance-min 1

~/medvedevGroup/dsk/build/bin/dsk2ascii -file ecoliReference.h5 -out ecoliReadsNew150-0-1ReferenceKmerCount.txt

echo "Kmer compare reference and bubblepopping"

python ~/medvedevGroup/bubblepopping/src/findCommonKmers.py ecoliReadsNew150-0-1BubblePoppedKmerCount.txt ecoliReadsNew150-0-1ReferenceKmerCount.txt > bubblepoppingConsoleOutputResultsKmerCompare.txt

echo "Kmer compare reference and bcalm"

python ~/medvedevGroup/bubblepopping/src/findCommonKmers.py ecoliReadsNew150-0-1BcalmKmerCount.txt ecoliReadsNew150-0-1ReferenceKmerCount.txt > bcalmConsoleOutputResultsKmerCompare.txt

