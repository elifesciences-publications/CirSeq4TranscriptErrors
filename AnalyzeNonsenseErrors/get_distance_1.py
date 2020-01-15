from __future__ import division
import sys
import numpy as np

#This script is to get,
#	1. the distance between the nonsense error and the original stop codon
#	2. the frequency of nonsense errors at the corresponding distance


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
			key2=Pos2
			if key2 in Pos_rate: 
				Pos_rate[key2].append(count_PTC[key]/cov_PTC[key])
			else:
				Pos_rate[key2]=[]
				Pos_rate[key2].append(count_PTC[key]/cov_PTC[key])

#Functions for the sliding-window analysis
def sw_analysis(rate_info):
	w=int(sys.argv[4])      #sliding_window size
	n=1000  #the upper length limit for bacterial mRNAs 
	count=[]
	tmp_list=[]
	for i in range(0,n+1):
		if i==0: #the 1st window
			for j in range(w):
				if j in rate_info:
					 count.append(sum(rate_info[j])/len(rate_info[j])) #Treat each gene equally
				else:
					count.append(0)
		else: #from the 2nd window
			if i-1 in rate_info:
				count.remove(sum(rate_info[i-1])/len(rate_info[i-1]))
			if i-1+w in rate_info:
				count.append(sum(rate_info[i-1+w])/len(rate_info[i-1+w]))
		tmp_list.append(sum(count)/w)
	return tmp_list

#Bootstrap
#Get a list for random sampling
list4rs=[]
for key in Pos_rate:
	for i in range(len(Pos_rate[key])): 
		list4rs.append(str(key)+'&&&'+str(i)) 

#Random sampling with replacement
n=len(list4rs)
tmp=np.random.choice(list4rs,(1000,n))
all_samples=tmp.tolist() # organize sampling resutls into one list

#Sliding-window analysis
all_results=[]	#store results for all sw_analysis
for i in all_samples:
	UPos_rate={}	#updatte Pos_rate for calculations 
	for j in i:
		info=j.split('&&&')
		key=int(info[0])
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
