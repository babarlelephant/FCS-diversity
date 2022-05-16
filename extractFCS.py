
g = open("fcs_AlphaCoV_Tegacovirus.tsv","w")
A = 16
B = A + 19

import os
os.system("c:\\code\\nextstrain_ncov\\blast\\bin\\makeblastdb -dbtype prot -in ref.fasta")
os.system("c:\\code\\nextstrain_ncov\\blast\\bin\\blastp -db ref.fasta -query sequence.fasta -evalue 1e-10 -outfmt 4 > a.txt")
f = open("a.txt","r")
def replaceDblBLanks(s):
	t = []
	b = False
	for i in range(len(s)):
		if not s[i] == " ":
			if not s[i] == "\r" and not s[i] == "\n":
				t.append(s[i])
			b = False
		else:
			if b == False:
				t.append(s[i])
				b = True
	return "".join(t).split(" ")

def extractFCS(ref,query):
	global A
	global B
	ref[1] = int(ref[1])
	ref[3] = int(ref[3])
	#print(ref[1],ref[3])
	if ref[1] <= A and ref[3] >= B:
		ind = ref[1]
		indA = -1
		indB = -1
		for i in range(len(ref[2])):
			if ind == A and indA == -1:
				indA = i
			if ind == B and indB == -1:
				indB = i
			if not ref[2][i] == "-":
				ind += 1
		return query[2][indA:indB+1]
	return ""
r=f.readline()

seqname = ""
query = ""
ref = ""

D = dict()

N = 0
while len(r) > 0:
	if r.startswith("Query= "):
		seqname= r.split("Query= ")[1]
	if r.startswith("Query_"):
		query = replaceDblBLanks(r)
	if r.startswith("0       "):
		ref = replaceDblBLanks(r)
		fcs = extractFCS(ref,query)
		if len(fcs) >0:
			N+=1
			if 0 and N == 1000:
				break
			if 0:
				g.write(">"+seqname)
				g.write(fcs+"\n")
				N+=1
				if N == 1000:
					quit()
			if 1:
				if not fcs in D:
					D.update({fcs:0})
				D[fcs]+=1
	r = f.readline()
	
t = []
for d in D:
	t.append([d,D[d]])

g.write("motif\tnb_sequences\n")

def fn(a):
	return -a[1]
t.sort(key=fn)
for a in t:
	g.write(a[0]+"\t"+str(a[1])+"\n")