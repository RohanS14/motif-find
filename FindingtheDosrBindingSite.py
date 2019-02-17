"""
DosR is a transcription factor that regulates the metabolic activity of Mycobacterium tubercolosis. In hypoxic conditions, it 
can persist in a latent state within the human body for many years. 
There are specific DNA sequences where DosR can bind, and the discovery of these binding sites has clinical significance. 
Here, the RandomizedMotifSearch algorithm has been used to find motifs that could be the DosR binding sites.

Randomized Algorithm for Motif Finding: a collection of randomly chosen k-mers from each string is chosen as the starting BestMotifs. 
By creating a Profile and finding the profile-most probable k-mer from each string(Motifs(Profile(motifs)etc.), 
we get a BestMotifs with lower and lower score, up to a point where it has found the local optimum. 
If RandomisedMotifSearch is run many times, we can choose the best one from many returned BestMotifs.
In particular, only if one of the implanted motifs is chosen as a random motif, then that may eventually lead the algorithm closer, 
or to the actual set of implanted motifs. If Dna was random then the profile would be completely uniform. 
But Dna is not random because of the higher probability of the motifs implanted,
and the profile points toward the implanted IF one of the implanted motifs is chosen by chance. That is its weakness, 
as it is guaranteed to reach the best scoring collection of Motifs that can be reached from its starting point (local optimum),
but if its randomly chosen starting point is wrong then it cannot reach the best solution from there, which is why it is run many times
to ensure that it finds a global, not just local optimum.

The results of the rpject are at the end.

"""

import random

def Count(Motifs):
	count = {}
	strLen = len(Motifs[0])
	for nucl in "ATCG":
		count[nucl] = []
		for colmn in range(strLen):
			count[nucl].append(0)
	strlistLen = len(Motifs)
	for strindex in range(strlistLen):
		for colmn in range(strLen):
			nucltd = Motifs[strindex][colmn]
			count[nucltd][colmn] += 1
	return count

def Profile(Motifs):
	strlistLen = len(Motifs)
	strLen = len(Motifs[0])
	count = Count(Motifs)
	profile = {}
	for nucl in "ATCG":
		profile[nucl] = []
		for colmn in range(strLen):
			profile[nucl].append(float(count[nucl][colmn])/strlistLen)
	return profile

def Consensus(Motifs):
	strLen = len(Motifs[0])
	#or Profile(Motifs)
	count = Count(Motifs)
	consensus = ""
	for colmn in range(strLen):
		m = 0
		fSymbol = ""
		for nucl in "ATCG":
			if count[nucl][colmn] > m:
				m = count[nucl][colmn]
				fSymbol = nucl
		consensus += fSymbol
	return consensus

def Score(Motifs):
	strLen = len(Motifs[0])
	strlistLen = len(Motifs)
	count = Count(Motifs)
	consensus = Consensus(Motifs)
	score = 0
	for colmn in range(strLen):
		for strindex in range(strlistLen):
			if Motifs[strindex][colmn] != consensus[colmn]:
				score += 1
	return score


def Pr(Text, Profile):
	p = 1#probability variable
	for nuclindex in range(len(Text)):
		p = p*Profile[Text[nuclindex]][nuclindex]
	return p

def ProfileMostProbablePattern(Text, k, Profile):
	probsDict = {}
	for nuclindex in range(len(Text)-k+1):
		probsDict[nuclindex] = [Pr(Text[nuclindex : nuclindex + k], Profile), Text[nuclindex : nuclindex + k]] 
	maxprob = -1
	maxkmer = ""
	for num in range(len(probsDict)):
		if probsDict[num][0] > maxprob:
			maxprob = probsDict[num][0]
			maxkmer = probsDict[num][1]
	return maxkmer



	BestMotifs = []
	for strindex in range(0, t):
		BestMotifs.append(Dna[strindex][0:k])
	strlen = len(Dna[0])
	for i in range(strlen-k+1):
		Motifs = []
		Motifs.append(Dna[0][i:i+k])
		for strindex in range(1, t):
			P = Profile(Motifs[0:strindex])
			Motifs.append(ProfileMostProbablePattern(Dna[strindex], k, P))
		if Score(Motifs) < Score(BestMotifs):
			BestMotifs = Motifs
	return BestMotifs

#Pseudocounts are used to avoid the issue of zero probability.
def CountWithPseudocounts(Motifs):
	count = {}
	strLen = len(Motifs[0])
	strlistLen = len(Motifs)
	for nucl in "ATCG":
		count[nucl] = []
		for colmn in range(strLen):
			count[nucl].append(1)
	for strindex in range(strlistLen):
		for colmn in range(strLen):
			nucltd = Motifs[strindex][colmn]
			count[nucltd][colmn] += 1
	return count

def ProfileWithPseudocounts(Motifs):
	strlistLen = len(Motifs)
	strLen = len(Motifs[0])
	count = CountWithPseudocounts(Motifs)
	profile = {}
	for nucl in "ATCG":
		profile[nucl] = []
		for colmn in range(strLen):
			profile[nucl].append(float(count[nucl][colmn])/(strlistLen+4))
	return profile

def Motifs(Profile, Dna):
	strLen = len(Dna[0])
	strlistLen = len(Dna)
	k = len(Profile["A"])
	Motifs = []
	for string in range(strlistLen):
		motifs = {}
		for i in range(strLen-k+1):
			motifs[i] = [Pr(Dna[string][i:i+k], Profile), Dna[string][i:i+k]]
			maxprob = -1
			maxkmer = ""
			for i in range(len(motifs)):
				if motifs[i][0] > maxprob:
					maxprob = motifs[i][0]
					maxkmer = motifs[i][1]
		Motifs.append(maxkmer)
	return Motifs

def RandomMotifs(Dna, k, t):
	M = len(Dna[0])-k
	RandMots = []
	for strindex in range(t):
		randpos = random.randint(1, (M+1))
		i = randpos - 1
		RandMots.append(Dna[strindex][i:i+k])
	return RandMots


def RandomizedMotifSearch(Dna, k, t):
	M = RandomMotifs(Dna, k, t)
	BestMotifs = M
	while True:
		Profile = ProfileWithPseudocounts(M)
		M = Motifs(Profile, Dna)
		if Score(M) < Score(BestMotifs):
			BestMotifs = M
		else:
			return BestMotifs


def RepeatedRMotifSearch(Dna, k, t, N):
	Dict = {}
	for r in range(N):
		Result =  RandomizedMotifSearch(Dna, k , t)
		Dict[r] = [Result, Score(Result)]
	minscore = 100000000
	minbestmots = []
	for i in range(len(Dict)):
		#print Dict[i][1], minscore 
		if Dict[i][1] < minscore:
			minscore = Dict[i][1]  
			minbestmots = Dict[i][0]
	return minbestmots

Dna = [   "GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

#DosR dormancy regulator binding sites 
t = 10 #no. of strings in DNA
k = 15 #length of binding site
N = 100 #no.of times algorithm is run
BestMotifs = RepeatedRMotifSearch(Dna, k, t, 100)
print BestMotifs 
print Score(BestMotifs)

"""Results of project

BestMotifs = ['GCGGACGAATGACCC', 'ATCGACCCGCGGCCC', 'ACCGTCGATGTGCCC', 'ATCGATCATCGGCCA', 'AAGGCCGAACGACCC',
	      'GAGGACCTTCGGCCC', 'GCGGACAAATGGCCC', 'TGGGACTTTCGGCCC','AAGGACTAACGGCCC', 'GGCCACCAATCGCCC']
Score(BestMotifs) = 43

"""
