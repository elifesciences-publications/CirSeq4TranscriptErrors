from __future__ import division
import sys
import numpy as np

#This script is to get,
#       1. the distance between the nonsense error and the original stop codon
#       2. the frequency of nonsense errors at the corresponding distance (further divided by the number of one-off codons)


#Extract info of the count and position of nonsense errors
count_PTC={}
file_P=open(sys.argv[1],'r')
for each_line in file_P:
	L=each_line.split()
 	key=L[0]+'&&&'+L[1]
	if key in count_PTC:
		count_PTC[key] +=1
	else:
		count_PTC[key]=1

#Extract info of the coverage of nonsense errors
cov_PTC={}
file_O=open(sys.argv[2],'r')
for each_line in file_O:
	L=each_line.split()
	ID=L[0]+'&&&'+L[1]
	if ID in count_PTC:
		if ID in cov_PTC:
			cov_PTC[ID] += int(L[4])
		else:
			cov_PTC[ID]=int(L[4])

#Get the distance to the original stop codon for nonsense errors & calculate the frequency
Pos_rate={}
Pos_ID={}
file_T=open(sys.argv[3],'r')
for each_line in file_T:
	L=each_line.split('\t')
	for key in count_PTC:
		info=key.split('&&&')
		if L[0]==info[0] and int(L[3])<= int(info[1]) <= int(L[4]):
			if L[6] in "+":
				Pos2=int(L[4])-int(info[1])-2
			elif L[6] in "-":
				Pos2=int(info[1])-int(L[3])-2
			key2=str(L[3])+'_'+str(L[4])	#cds_ID
			cds_length=int(L[4])-int(L[3])+1
			rate=count_PTC[key]/cov_PTC[key]
			value=str(cds_length)+'&&&'+str(Pos2)+'&&&'+str(rate)	#length_dis_rate
			if key2 in Pos_rate:
				Pos_rate[key2].append(value)
			else:
				Pos_rate[key2]=[]
				Pos_rate[key2].append(value) 

#Get the sequence info from cds.fa
cds_fa={}
seq=""
file_P=open(sys.argv[4],'r')
for each_line in file_P:
        if each_line[0] == ">":
                if seq != "":
                        cds_fa[ID]=seq
                        seq=""
		start=int((each_line.split(':')[1].strip()).split('-')[0])+1
		end=str((each_line.split('-')[1].strip()).split('(')[0])
                ID=str(start)+'_'+str(end)
        else:
                seq += each_line.strip()
        if seq != "":
                cds_fa[ID]=seq
#Define one-off codons
one_offs=['AAA','CAA','GAA','TTA','TCA','TGA','TAT','TAC','TAG','AAG','CAG','GAG','TTG','TCG','TGG','TAA','TAT','TAC','AGA','CGA','GGA','TAA','TTA','TCA','TGT','TGC','TGG']

#Define the function to get # of oneoffs
def oneoff(pos,cds_id,w_s):
	if pos%3 ==0:	#5' extension
		focus_seq=cds_fa[cds_id][-(w_s+pos+2):-(pos+1)]+cds_fa[cds_id][-(pos+1)]	
	elif pos%3 ==1:	# 5' and 3' extesnsions
		pos -= 1
		focus_seq=cds_fa[cds_id][-(w_s+pos+2):-(pos+1)]+cds_fa[cds_id][-(pos+1)]
	elif pos%3 ==2:	# 3' extension
		pos -= 2
		focus_seq=cds_fa[cds_id][-(w_s+pos+2):-(pos+1)]+cds_fa[cds_id][-(pos+1)]
	Num_codon=len(focus_seq)/3
	Num_oneoff=0
	for h in range(int(Num_codon)):
		if focus_seq[h*3:h*3+2]+focus_seq[h*3+2] in one_offs:
			Num_oneoff +=1
	return Num_oneoff

#Functions for sliding window analysis
def sw_analysis(rate_info):
	w=int(sys.argv[5])      #sliding_window size
	n=1000
	tmp_list=[]
	for i in range(0,n+1):
		count=[]
		for j in range(i,i+w): #till i+w-1
			tmp1=[]
			for key in rate_info:
				for m in rate_info[key]:
					info1=m.split('&&&')
					if int(info1[1]) == j:
						Num_oneoff=oneoff(i,key,w)
						tmp1.append(float(info1[2])/Num_oneoff)
			if len(tmp1) >0:
				count.append(sum(tmp1)/len(tmp1))
			else:
				count.append(0)
		tmp_list.append(sum(count)/w)       #Treat each site equally
	return tmp_list 

#Bootstrap
#Get a list for random sampling
list4rs=[]
for key in Pos_rate:
	for i in range(len(Pos_rate[key])): 
		list4rs.append(str(key)+'***'+str(i)) 

#Random sampling with replacement
n=len(list4rs)
tmp=np.random.choice(list4rs,(1000,n))
all_samples=tmp.tolist() # organize sampling resutls into one list

#Calculate for each sample
all_results=[]	#store results for all sw_analysis
for i in all_samples:
	UPos_rate={}	#updatte Pos_rate for calculations 
	for j in i:
		info=j.split('***')
		key=info[0]
		value=Pos_rate[key][int(info[1])]
		if key in UPos_rate:
			UPos_rate[key].append(value)
		else:
			UPos_rate[key]=[]
			UPos_rate[key].append(value)
	all_results.append(sw_analysis(UPos_rate)) #Ultimately get a 1000 X len(list4rs) maxtrix

#Output
print 'Distance'+'\t'+'Frequency'+'\t'+'StandardDeviation'+'\t'+'Sample'
out_results=sw_analysis(Pos_rate)	# the same length as each element in all_results
for i in range(len(out_results)):
	pos=i
	value=out_results[i]
	tmp4std=[]	
	for j in range(len(all_results)):
		tmp4std.append(all_results[j][i])
	std=np.std(tmp4std)
	print str(i)+'\t'+str(value)+'\t'+str(std)+'\t'+"Samples"
