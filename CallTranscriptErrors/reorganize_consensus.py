import sys,re

outfile1=open(sys.argv[1]+"_"+"reads.fq",'w')	#Reorganized reads

#Get the 1st hit of each read
D1={}
file_1=open(sys.argv[1],'r')
for each_line1 in file_1:
	if each_line1[0] !='@' and each_line1.split()[0] not in D1:
		key=each_line1.split()[0]
		value=each_line1.split()[1:]
		D1[key]=value

#Reconstruct the consensus sequence
def parse_cigar(cigar):
	pattern = re.compile('([MIDNSHPX=])')
	values = pattern.split(cigar)[:-1]
	paired = (values[n:n+2] for n in xrange(0, len(values), 2)) 
	M,S,I=[],[],[0]
	for pair in paired:
		t=pair[1] 
		l=int(pair[0])  
		if t == 'M':
			M.append(l)
		elif t== 'S':
			S.append(l)
		elif t=='I':
			I.append(l)
	return(M,S,I)

for key in D1:   #Analyze Cigar information
	cigar=D1[key][4]
	(M,S,I)=parse_cigar(cigar)  # M can have only one value and S can have two.
	if len(M)*len(S)>0: # make sure to have mapping information 
		l=len(D1[key][8])/2
		mismatch=sum(S)
		match=sum(M)
		insertion=sum(I) 
		segment=l-M[0]
		if abs(match-mismatch) <=8: 
			if len(M)==len(S)==1:
				seq=D1[key][8][:l]
				qs=D1[key][9][:l]
			else:
				seq=D1[key][8][S[0]:(S[0]+l)]
				qs=D1[key][9][S[0]:(S[0]+l)]
			outfile1.write("@"+key+"\n"+seq+"\n"+"+"+"\n"+qs+"\n")
		elif 24>=(mismatch-match)>8:
			if len(S)==2: 
				seq=D1[key][8][S[0]:(S[0]+match+insertion)]
				qs=D1[key][9][S[0]:(S[0]+match+insertion)]
				adjust=l-(match+insertion)
				if (S[0]-adjust) >=0: 
					barcode=D1[key][8][(S[0]-adjust):S[0]]  
				else:
					barcode=D1[key][8][S[0]-adjust]+D1[key][8][:S[0]]
			elif len(S)==1:
				seq=D1[key][8][:(match+insertion)]
				qs=D1[key][9][:(match+insertion)]
				barcode=D1[key][8][(match+insertion):l]
			outfile1.write("@"+key+"\n"+seq+"\n"+"+"+"\n"+qs+"\n")   
