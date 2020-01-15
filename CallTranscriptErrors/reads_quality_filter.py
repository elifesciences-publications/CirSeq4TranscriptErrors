import sys

# mark low quality calls as "N" and filter out supplementary aligments

#List of base call qualities
delete_list=['"','#','$','%','&','\'','(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I','J','K']


D1={}  #This is used to filter out supplementary alignments
file_P=open(sys.argv[1],'r')
for each_line in file_P:
	if each_line[0] =='@':
		print each_line.strip()
	else:
		if each_line.split()[5] !="*" and each_line.split()[0] not in D1:#focus on reads with hits
			D1[each_line.split()[0]]=''	#Just output the 1st(best) read alignment
			low_qual_list=[]
			L=each_line.split()
			qual=list(L[10])
			reads=list(L[9])
			for i in range(len(qual)):
				if qual[i] in delete_list: #mark qs that are < 70. Info on self dedined qs scale can be found in get_consensus.py
					low_qual_list.append(i)
				else:
					qual[i]='K'
			for i in low_qual_list:
				reads[i]='N'
			L[9]=''.join(reads)
			L[10]=''.join(qual)
			print '\t'.join(L)
