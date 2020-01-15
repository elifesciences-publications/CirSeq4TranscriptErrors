from __future__ import division
import sys

#read the fastq file
def parse_fastq(filename):
	file = open(filename)
	R={} #Extract reads
	for i,line in enumerate(file):
		if i % 4 == 0:
			key=line.split()[0]
		elif i % 4 == 1:
			R[key]=line.rstrip('\n')
	return R

R1=parse_fastq(sys.argv[1])

#Select reads longer than 60
list=[]
for key in R1:
	r1=R1[key]
	if len(r1) <60:
		list.append(key)
for i in list:
	del R1[i]

#Infer the size of tandem repeats
cutoff=float(sys.argv[2]) #cutoff on the identify of repeats 
for key in R1:
	r1=R1[key]
	L=len(r1)    
	C={}
	for p in range(10,L-49):    # L-p is the length of sites that are checked --> repeats with length larger than 250 can not be detected
		c1=0
		for i in range(0,L-p):	#this allows to evaluate more than one repeat 
			if r1[i]==r1[i+p]:
				c1+=1
		C[p]=c1/(L-p)
	highest=max(C.values())
	if highest >=cutoff:
		P_list=[k for k,v in C.items() if v == highest]
		print key,"\t",min(P_list)
