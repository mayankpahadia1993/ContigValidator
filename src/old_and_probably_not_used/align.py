'''

	Author - Mayank Pahadia
	Program - Align the fasta file to the reference genome

	Run the program - python align.py referenceGenome.fa inputFasta.fa outputResults.txt


'''

from suffix_tree import SuffixTree
import sys

def revc(a):
	b=a[::-1]
	c=''
	for i in range(0,len(b)):
		if(b[i]=='A' or b[i]=='a'):
			c+='T'
		elif(b[i]=='C' or b[i]=='c'):
			c+='G'
		elif(b[i]=='G' or b[i]=='g'):
			c+='C'
		elif(b[i]=='T' or b[i]=='t'):
			c+='A'
	return c

def find_substring(stree, substring):
	for i in stree.leaves:
		j = i.pathLabel.replace('$','')
		if j==substring:
			return 1
	return 0

if(len(sys.argv)<4):
	print "Please pass the reference genome, inputfile and outputfile to the python script"
	print "Example - python align.py referenceGenome.fa inputFasta.fa outputResults.txt"
	exit()

referenceGenomeFile = sys.argv[1]
inputFastaFile = sys.argv[2]
outputResultsFile = sys.argv[3]
referenceGenome=[]
inputFasta=[]
inputFastaId=[]
# outputResults={}
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
with open(inputFastaFile) as f:
	temp=""
	for line in f:
		line=line.replace('\n','')
		line = line.upper()
		if(line==""):
			continue
		if(line[0]=='>'):

			if(temp!=""):
				inputFasta.append(temp)
			temp=""
			# outputResults[line]=0
			inputFastaId.append(line)
			continue
		temp+=line
inputFasta.append(temp)
print "Input Fasta read"
# print "inputFasta coming up"
# for i in inputFasta:
# 	print i


if(len(inputFasta)!=len(inputFastaId)):
	print "some problem because total number of reads don't match total number of ids"
	exit()
print "Suffix Tree start"
stree = SuffixTree(referenceGenome[0])
print "Suffix Tree created"
with open(outputResultsFile,'w') as f:
	for i in range(0,len(inputFasta)):
		# outputResults[inputFastaId[i]] = find_substring(stree,inputFasta[i])
		# c = inputFastaId[i] + d + str(find_substring(stree,inputFasta[i]))+ "\n"
		if(stree.find_substring(inputFasta[i]) >= 0 or stree.find_substring(revc(inputFasta[i])) >= 0):
			find=1
			countAlignments+=1
		else:
			find=0
		c = inputFastaId[i] + d + str(find)+ "\n"
		f.write(c)

print "total Alignments - ",countAlignments
percent = 1.0*countAlignments/len(inputFastaId)
percent = percent*100.0
print "percentage of Alignment - ", percent  

exit()