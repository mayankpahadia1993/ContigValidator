echo "Alignment of Reference to bcalm output"

python ~/medvedevGroup/bubblepopping/src/align.py ../ecoliReference.fa ecoliReadsNew150-0-1.unitigs.fa bcalmResultsSuffix.txt > bcalmConsoleOutputResultsAlignment.txt

echo "Alignment of Reference to bubblepopping output"

python ~/medvedevGroup/bubblepopping/src/align.py ../ecoliReference.fa ecoliReadsNew150-0-1BubblePoppedFastaFile.fa bubblePoppedResultsSuffix.txt > bubblepoppingConsoleOutputResultsAlignment.txt


