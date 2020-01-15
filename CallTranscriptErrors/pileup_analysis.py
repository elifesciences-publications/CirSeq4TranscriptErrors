from __future__ import division
import sys,re

L=['+','-']
indel_count=0
cutoff=float(sys.argv[2])	#cutoff of the error frequency

#header 
print 'Chr/Gene','\t','Pos','\t','Ref','\t','Sub','\t','Total_coverage','\t','Count_of_errors'

file_P=open(sys.argv[1],'r')
for each_line in file_P:
	if int(each_line.split()[3])>0:
		comma,un,nun=0,0,0
		count_dic={'A':0,'T':0,'C':0,'G':0,'a':0,'t':0,'c':0,'g':0}
		sum_dic={}
		bases=each_line.split()[4]
		total=int(each_line.split()[3])
		if any(x in bases for x in L):
			indel_count+=1
		else:	#detect only substitutions 
			for i in range(len(bases)): 
				if bases[i] == ',':
					comma+=1
				elif bases[i] in ('A','T','C','G','a','t','c','g'):
					count_dic[bases[i]]+=1
				elif bases[i] in ('N','n'):
					un+=1
				elif bases[i] == '^':	# ^ marks the start of a read segment and is followed by quality base; $ marks the end of a read segment and no quality base follows. 
						if bases[i+1] in ('A','T','C','G','a','t','c','g'):
							count_dic[bases[i+1]]-=1
						if bases[i+1] in ('N','n'):
							nun+=1
			Un=un-nun
			total-=Un	# exclude Ns and get the effective coverage
			if total >0: 
				Mut=sum(count_dic.itervalues())
				sum_dic['A']=count_dic['A']+count_dic['a'] 
				sum_dic['T']=count_dic['T']+count_dic['t']
				sum_dic['C']=count_dic['C']+count_dic['c']
				sum_dic['G']=count_dic['G']+count_dic['g']
				num_types=0
				mut={}
				for key in sum_dic:
					if sum_dic[key]>0: 
						num_types+=1
						mut[key]=sum_dic[key]
				if 0< Mut/total<=cutoff:
					for key in mut:
						for i in range(int(mut[key])):  
							print each_line.split()[0],'\t',each_line.split()[1],'\t',each_line.split()[2],'\t',key,'\t',total,'\t',mut[key] 
