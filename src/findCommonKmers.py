'''

	Author - Mayank Pahadia
	Program - Find the common kmers between two kmercount files. 

	Run the program - python findCommonKmers.py input1.txt input2.txt


'''

import sys
from collections import defaultdict

if(len(sys.argv)<5):
	print "Please pass the input files"
	print "Example - python findCommonKmers.py tempKmerOut.txt input1.txt input2.txt commonKmersFrom1in2.txt commonKmersFrom2in1.txt"
	exit()

tempKmerOut=sys.argv[1]
input1 = sys.argv[2]
inputs = sys.argv[3:]

kmers = defaultdict(lambda: 0)

print input1, inputs


countInput1=0

delim="\t"


out = open(tempKmerOut,'w')
c=""
c="Percentage of Kmers from Reference File in Input File" + delim + "Percentage of Kmers from Input File in Reference File" + "\n"
out.write(c)

with open(input1) as f:
	for line in f:
		line=line.replace('\n','')
		line=line.split()
		kmers[line[0]]=1
		countInput1+=1

c=''
ijk=0
print countInput1
while(ijk<len(inputs)):
	commonFrom1 = defaultdict(lambda: 0)
	count=0
	countInput2=0
	input2 = inputs[ijk]
	ijk+=1
	file12=open(inputs[ijk],'w')
	ijk+=1
	file21=open(inputs[ijk],'w')
	ijk+=1
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

	# print "total kmers in", input1, " - ", countInput1
	# print "total kmers in", input2, " - ", countInput2
	# print "common kmers - ", count
	if(countInput1 == 0):
		percent1=0
	else:
		percent1 = count*1.0 / countInput1
		percent1*=100.0

	# print "common kmers from the total kmers in", input1, " - ", percent, "%"
	
	if(countInput2 == 0):
		percent2=0
	else:
		percent2 = count*1.0 / countInput2
		percent2*=100.0
	

	# print "common kmers from the total kmers in", input2, " - ", percent, "%"

	c=""
	c+=str(percent1) + "%" + delim + str(percent2) + "%" + "\n"
	out.write(c)

	file12.close()
	file21.close()

out.close()
