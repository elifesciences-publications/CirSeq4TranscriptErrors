import sys

#Get ID of reads with MAPQ==0(no/multiple hits)
ID_list={}
file_1=open(sys.argv[1],'r')
for each_line1 in file_1:
	if each_line1[0] !='@' and each_line1.split()[4] == '0':
		ID_list["@"+each_line1.split()[0]]=""

#Filter out above IDs in the fasta file
file_2=open(sys.argv[2],'r')  
for each_line2 in file_2: 
	if each_line2.rstrip('\n') in ID_list:
			file_2.next()
			file_2.next()
			file_2.next()
	else: 
		print each_line2.rstrip('\n')
