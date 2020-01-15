from __future__ import division
import sys,operator,math

outfile1=open(sys.argv[1]+'_'+"Out.check",'w') #Output for double check: ID+sequence alignment+ QS alighment+Consensus_seq+Consensus_qs
outfile2=open(sys.argv[1]+'_'+"Consensus.fq",'w')   #Consensus sequence (duplicates) in FASTQ format

def parse_fastq(filename):
	file = open(filename)
	R,Q={},{} #Extract reads and quality scores
	for i,line in enumerate(file):
		if i % 4 == 0:
			key=line.split()[0]
		elif i % 4 == 1:
			R[key]=line.rstrip('\n')
		elif i % 4 == 3:
			Q[key]=line.rstrip('\n') 
	return (R,Q)

def read_i(filename): # Read info on the size of tandem repeats
	I={}
	file_I=open(filename,'r')
	for each_line in file_I:
		key=each_line.split()[0]
		value=int(each_line.split()[1])
		I[key]=value
	return I

def align_matrix(D1,P): #Construct the 'matrix' to construct consensus duplicates
	infor_matrix={} 
	for key in D1:
		p=P[key]
		L2=len(D1[key])
		n2=L2//p
		m2=L2%p
		infor_list=[]
		for i in range(0,n2):
			seq=D1[key][i*p:(i+1)*p]
			infor_list.append(seq)
		if m2>0:
			seq=D1[key][n2*p:]+(p-m2)*'N'
			infor_list.append(seq)
		infor_matrix[key]=infor_list  
	return infor_matrix

#Read fastq
(R1,Q1)=parse_fastq(sys.argv[1])

#Read info on the size of tandem repeats 
P=read_i(sys.argv[2])

#Parse reads with infered tandem repeats
r_list=[]
for key in R1:
	if key not in P:
		r_list.append(key)
for i in r_list:
	del R1[i],Q1[i]

#Construct the 'matrix' to construct consensus duplicates
seq_matrix=align_matrix(R1,P)
quality_matrix=align_matrix(Q1,P)

#Transfer ASCII(Illumina v1.8 and later) to quality scores
transfer_dic={'"':1,'#':2,'$':3,'%':4,'&':5,'\'':6,'(':7,')':8,'*':9,'+':10,',':11,'-':12,'.':13,'/':14,'0':15,'1':16,'2':17,'3':18,'4':19,'5':20,'6':21,'7':22,'8':23,'9':24,':':25,';':26,'<':27,'=':28,'>':29,'?':30,'@':31,'A':32,'B':33,'C':34,'D':35,'E':36,'F':37,'G':38,'H':39,'I':40,'J':41,'K':42}

transfer_back_dic=dict((v,k) for k,v in transfer_dic.iteritems())

#Get the consensus
Con={}
re_qs_matrix={}
for key in seq_matrix:
	KEY=key
	p=P[key]
	read_list=seq_matrix[key]
	Consensus_seq=''
	for i in range(0,p):
		base_dic={'A':0, 'T':0, 'C':0, 'G':0}
		for j in range(len(read_list)):
			if read_list[j][i] in ('A','T','C','G'):
				base_dic[read_list[j][i]]+=1
		max_base=max(base_dic.iteritems(), key=operator.itemgetter(1))[0] # identify the base(s) with highest counts. If there are two maxium values, will be filter out in the next step.
		total=base_dic['A']+base_dic['T']+base_dic['C']+base_dic['G']
		if total >=2:	#At least 2 infromative bases covered 
			if base_dic[max_base]/total > 0.65 and base_dic[max_base] >=2: # 0.65 is arbitrary 
				Consensus_seq+=max_base
			else:
				Consensus_seq+='N' 
		else:
			Consensus_seq+='N'
	Con[key]=Consensus_seq
#recalculate QS
	qs_list=quality_matrix[key]
	re_qs_list=[]
	for i in range(0,p):
		if Consensus_seq[i]=='N':
			re_qs_list.append('"')  # assign the lowest qs
		else:
			qs_dic={'A':1, 'T':1, 'C':1, 'G':1}  #reacord probabilities
			for j in range(len(qs_list)):
				for h in ('A','T','C','G'):
					if read_list[j][i] in ('A','T','C','G'):
						if read_list[j][i]==h:
							qs_dic[h]*=1-10**(transfer_dic[qs_list[j][i]]/(-10))
						else:
							qs_dic[h]*=(10**(transfer_dic[qs_list[j][i]]/(-10)))/3
			denominator=0
			for key in qs_dic:
				if key ==Consensus_seq[i]:
					numerator=qs_dic[key]
					denominator+=qs_dic[key]  
				else:
					denominator+=qs_dic[key]
			re_qs_list.append(numerator/denominator)
	re_qs_matrix[KEY]=re_qs_list

#Output in FASTQ format(reform qs lines)
for key in re_qs_matrix:
	for i in range(len(re_qs_matrix[key])):
		if re_qs_matrix[key][i]=='"':
			re_qs_matrix[key][i]=='"'
		elif re_qs_matrix[key][i] == 1:
			re_qs_matrix[key][i] = 'K'    #assign highest qs
		else:
			Phred=round(-10*math.log10(1-re_qs_matrix[key][i]))  #1-proability of a correct call
			if Phred >= 70:
				re_qs_matrix[key][i] = 'W'	#a self-defined notation that corresponds to a quality score >= 70
			elif 42< Phred<70:
				re_qs_matrix[key][i] = 'K'
			elif 0< Phred<=42:  
				tem=int(Phred)
				re_qs_matrix[key][i]=transfer_back_dic[tem]
			elif Phred==0:
				re_qs_matrix[key][i]='"'
	outfile2.write(key+'\n'+2*Con[key]+'\n'+'+'+'\n'+2*(''.join(re_qs_matrix[key]))+'\n')

#Output multiple alignment + consensus sequence for doublechecks and adjusting parameters
for key in seq_matrix:
	outfile1.write(key+'\n')
	for i in seq_matrix[key]:
		outfile1.write(i+'\n')
	for i in quality_matrix[key]:
		outfile1.write(i+'\n')
	outfile1.write(Con[key]+'\n'+''.join(re_qs_matrix[key])+2*'\n')
