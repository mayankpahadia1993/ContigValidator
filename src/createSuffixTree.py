'''

	Author - Mayank Pahadia
	Program - Create a suffix tree and store it in a file

	Run the program - python createSuffixTree.py referenceGenome.fa suffixTreeOutput.p


'''

from suffix_tree import SuffixTree
import sys
try:
    import cPickle as pickle
except:
    import pickle 

if(len(sys.argv)<3):
	print "Please pass the reference genome, and outputfile to the python script"
	print "Example - python createSuffixTree.py referenceGenome.fa suffixTreeOutput.p"
	exit()

referenceGenomeFile = sys.argv[1]
outputResultsFile = sys.argv[2]
referenceGenome=[]

temp=""
d="\t"
countAlignments=0
with open(referenceGenomeFile) as f:
	temp=""
	for line in f:
		line=line.replace('\n','')
		line = line.upper()
		if(line==""):
			continue
		if(line[0]=='>'):
			if(temp!=""):
				referenceGenome.append(temp)
			temp=""
			continue
		temp+=line
referenceGenome.append(temp)
# print referenceGenome[0]
print "Reference Genome read"


print "Suffix Tree start"
stree = SuffixTree(referenceGenome[0])
print "Suffix Tree created"

pickle.dump( stree, open( outputResultsFile, "wb" ) )

print "Suffix Tree saved"

exit()

