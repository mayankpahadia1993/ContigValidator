RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' #No Color
cd ..
bash run.sh -r test/reference.fa -r test/reference2.fa -s test/suffixtree.p -i test/reads.fa -a test/alignresults.txt -kmer-size 25 -abundance-min 1 2> tempError.txt
if [ "$?" = 0 ]; then
	printf "${BLUE}TEST PASSED SUCCESSFULLY${NC}\n"
else
	printf "${RED}ERROR - "
	# echo -e "I ${RED}love${NC}"
	cat tempError.txt

	printf "${NC}\n"
fi

rm -f tempError.txt

# find test/ -type f ! -name 'reference.fa' ! -name 'reads.fa' ! -name 'test.sh' ! -name 'alignresults.txt' -delete