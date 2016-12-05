cd ..
bash run.sh -r test/reference.fa -s test/suffixtree.p -i test/reads.fa -a test/alignresults.txt -abundance-min 1
find test/ -type f ! -name 'reference.fa' ! -name 'reads.fa' ! -name 'test.sh' ! -name 'alignresults.txt' -delete