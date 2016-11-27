'''

	Author - Mayank Pahadia
	Program - Find the common kmers between two kmercount files. 

	Run the program - python findCommonKmers.py input1.txt input2.txt


'''

import sys
from collections import defaultdict

if(len(sys.argv)<5):
	print "Please pass the input files"
	print "Example - python findCommonKmers.py input1.txt input2.txt commonKmersFrom1in2.txt commonKmersFrom2in1.txt"
	exit()

input1 = sys.argv[1]
input2 = sys.argv[2]
file12=open(sys.argv[3],'w')
file21=open(sys.argv[4],'w')
kmers = defaultdict(lambda: 0)
commonFrom1 = defaultdict(lambda: 0)
count=0
countInput1=0
countInput2=0
delim="\t"
with open(input1) as f:
	for line in f:
		line=line.replace('\n','')
		line=line.split()
		kmers[line[0]]=1
		countInput1+=1

c=''
with open(input2) as f:
	for line in f:
		line=line.replace('\n','')
		countInput2+=1
		line=line.split()
		c=''
		c+=line[0] + delim + str(kmers[line[0]])
		c+='\n'
		file12.write(c)
		if(kmers[line[0]]==1):
			commonFrom1[line[0]]=1
			count+=1

for key in kmers:
	c=''
	c+=str(key) + delim
	if commonFrom1[key]==1:
		c+=str(commonFrom1[key])
	else:
		c+=str(0)
	c+='\n'
	file21.write(c)

print "total kmers in", input1, " - ", countInput1
print "total kmers in", input2, " - ", countInput2
print "common kmers - ", count

percent = count*1.0 / countInput1
percent*=100.0

print "common kmers from the total kmers in", input1, " - ", percent, "%"

percent = count*1.0 / countInput2
percent*=100.0

print "common kmers from the total kmers in", input2, " - ", percent, "%"

file12.close()
file21.close()


