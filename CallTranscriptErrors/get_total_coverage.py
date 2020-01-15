from __future__ import division
import sys

freq=float(sys.argv[2])
cutoff=int(1/freq)

file_R=open(sys.argv[1],'r')
S,indel_count=0,0
INDEL=['+','-']
Count={}
for each_line in file_R:
	if int(each_line.split()[3])>0:
		seq=each_line.split()[4]
		if any(x in seq for x in INDEL):
			indel_count+=1
		else:
			total=int(each_line.split()[3])
			seq=each_line.split()[4]
			base=each_line.split()[2]
			un,nun=0,0
			count_dic={'A':0,'T':0,'C':0,'G':0,'a':0,'t':0,'c':0,'g':0}
			sum_dic={}
			if total >=cutoff:
				for i in range(len(seq)):
					if seq[i] in ('A','T','C','G','a','t','c','g'):
						count_dic[seq[i]]+=1
					elif seq[i] in ('N','n'):
						un+=1
					elif seq[i] == '^':
						if seq[i+1] in ('A','T','C','G','a','t','c','g'):
							count_dic[seq[i+1]]-=1
						if seq[i+1] in ('N','n'):
							nun+=1
				Un=un-nun
				total-=Un	# exclude Ns and get the effective coverage    
				if total >=cutoff:
					sum_dic['A']=count_dic['A']+count_dic['a']
					sum_dic['T']=count_dic['T']+count_dic['t']
					sum_dic['C']=count_dic['C']+count_dic['c']
					sum_dic['G']=count_dic['G']+count_dic['g']
					#num_types=0
					#mut={}
					#for key in sum_dic:
					#	if sum_dic[key]>0: 
					#		num_types+=1
				if total >=cutoff:
					S+=total	# total coverages
					if base in Count:	# for coverages of four bases 
						Count[base]=Count[base]+total
					else:
						Count[base]=total
print "Total Coverage:",S  

for i in Count:
	print i,"\t",Count[i]
