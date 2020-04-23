import sys,re
#consider hard clip scenarios 
#consider the length of one repeat

outfile1=open(sys.argv[1]+"_"+"reads_indel.fq",'w')


#Pick up the 1st hits of each duplicate reads
D1={}
file_1=open(sys.argv[1],'r')
for each_line1 in file_1:
	if each_line1[0] !='@' and each_line1.split()[0] not in D1:
		key=each_line1.split()[0]
		value=each_line1.split()[1:]
		D1[key]=value
#print D1 

#Resolve the read end points & barcodes

def parse_cigar(cigar):
	pattern = re.compile('([MIDNSHPX=])')
	values = pattern.split(cigar)[:-1]
	paired = (values[n:n+2] for n in xrange(0, len(values), 2))
	M,S,I,D,H=[],[],[],[],[]
	for pair in paired:
		t=pair[1]
		l=int(pair[0])
		if t == 'M':
			M.append(l)
		elif t== 'S':
			S.append(l)
		elif t=='I':
			I.append(l)
		elif t=='D':
			D.append(l)
		elif t=='H':
			H.append(l)
	return(M,S,I,D,H)

#Extract mapping pattern from cigar
def extract_cigar(cigar):
	pattern=''
	for i in cigar:
		if i =="S" or i=="M" or i=="H":
			pattern+=i
	return(pattern)

#Go through each reads
for key in D1:   #Analyze Cigar information
	cigar=D1[key][4]
	(M,S,I,D,H)=parse_cigar(cigar)  #
	pattern=extract_cigar(cigar)
	L=len(D1[key][8])/2  
	if len(M)==2 and len(I)==1 and len(D)==0:  # For potential insertions
		l=M[0]+M[1]+I[0]	#
		if L-5 <= l <= L+0:	# 5: arbitary 
			if pattern =='MMS':
				seq=D1[key][8][:l]
				qs=D1[key][9][:l]
			elif pattern =='SMMS' or pattern =='SMM': # 'SMM' may not be regular repeats but try to recover as much reads as possible.
				seq=D1[key][8][S[0]:(S[0]+l)]
				qs=D1[key][9][S[0]:(S[0]+l)]
			elif pattern =='MMH' or pattern =='HMMH' or pattern =='HMM':	#potential insertions with hard clips
				seq=D1[key][8][:l]
				qs=D1[key][9][:l]
			outfile1.write("@"+key+"\n"+seq+"\n"+"+"+"\n"+qs+"\n")
	elif len(M)==2 and len(D)==1 and len(I)==0:  # For potential deletions
		l=M[0]+M[1]	#exclude D
		if L-5 <= l <= L+0:
			if pattern =='MMS':
				seq=D1[key][8][:l]
				qs=D1[key][9][:l]
			elif pattern =='SMMS' or pattern =='SMM':
				seq=D1[key][8][S[0]:(S[0]+l)]
				qs=D1[key][9][S[0]:(S[0]+l)]
			elif pattern =='MMH'or pattern =='HMMH' or pattern =='HMM':
				seq=D1[key][8][:l]
				qs=D1[key][9][:l]
			outfile1.write("@"+key+"\n"+seq+"\n"+"+"+"\n"+qs+"\n")
