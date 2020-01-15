import sys

def parse_fastq(filename):
	file = open(filename)
	R={} #Extract reads
	for i,line in enumerate(file):
		if i % 4 == 0 or i % 4 == 2:
			print line.rstrip('\n')
		elif i % 4 == 1 or i % 4 == 3:
			L=line.rstrip('\n')
			if len(L)>8:
				print L[4:-4]
parse_fastq(sys.argv[1])
