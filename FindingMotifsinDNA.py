import random

#Three different algorithms to find regulatory motifs in DNA.

#Returns a Count matrix for each nucleotide in each column
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

""" Example:
	  Motifs = ["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC", "TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT", "GGCAAGTACC"]
		Count(Motifs) = {'A': [2, 3, 3, 3, 6, 4, 2, 2, 1, 3], 
		 'C': [2, 3, 4, 3, 2, 3, 2, 1, 3, 3], 
		 'T': [2, 2, 0, 4, 1, 0, 2, 2, 1, 4], 
		 'G': [4, 2, 3, 0, 1, 3, 4, 5, 5, 0]}"""


#The function Profile(motifs) is very similar to Count(motifs), except sum of values in each column adds up to 1(as a probability)
#Input:a set of slightly varying motifs
#Output: the probability of the occurrence of each nucleotide at each position of the strings based on observed occurrences
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

"""Example:
	print Profile(["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC", "TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT", "GGCAAGTACC"])
	{'A': [0.2, 0.3, 0.3, 0.3, 0.6, 0.4, 0.2, 0.2, 0.1, 0.3], 
 	'C': [0.2, 0.3, 0.4, 0.3, 0.2, 0.3, 0.2, 0.1, 0.3, 0.3], 
 	'T': [0.2, 0.2, 0.0, 0.4, 0.1, 0.0, 0.2, 0.2, 0.1, 0.4], 
 	'G': [0.4, 0.2, 0.3, 0.0, 0.1, 0.3, 0.4, 0.5, 0.5, 0.0]}"""

#The consensus motif, "correct" motif
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

"""Example:
print Consensus(["GTACAACTGT", 
				 "CAACTATGAA", 
				 "TCCTACAGGA", 
				 "AAGCAAGGGT", 
				 "GCGTACGACC", 
				 "TCGTCAGCGT", 
				 "AACAAGGTCA", 
				 "CTCAGGCGTC", 
				 "GGATCCAGGT", 
				 "GGCAAGTACC"])
    Output: GACTAAGGGT
"""

#Computes the total number of "wrong" nucleotides in each column (vs. consensus string)
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
	
# Example:- Score(["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA",  "AAGCAAGGGT",  "GCGTACGACC",  "TCGTCAGCGT",  "AACAAGGTCA",  "CTCAGGCGTC",  "GGATCCAGGT",  "GGCAAGTACC"]) = 57

#Using a profile, calculates the probability of a motif occurring
def Pr(Text, Profile):
	p = 1#probability variable
	for nuclindex in range(len(Text)):
		p = p*Profile[Text[nuclindex]][nuclindex]
	return p

"""Example:
	Prof = {"A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9], "C" : [0.2, 0.3, 0.0, 0.4, 0.0, 0.1], "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0], "T" : [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]}
	Pr("TCGGTA", Prof) = 0.00405"""

#Returns the first occurring most similar motif to a certain string of length k (k-mer) according to probability from the k-mer's profile  
def ProfileMostProbablePattern(Text, k, Profile):
	probsDict = {}
	for nuclindex in range(len(Text)-k+1):
		probsDict[nuclindex] = [Pr(Text[nuclindex : nuclindex + k], Profile), Text[nuclindex : nuclindex + k]] 
		#In the dictionary ProbsDict: key is k to maintain reference to numerical order, value is a list with 0 item as probability and 1 item as k-mer
	maxprob = -1
	maxkmer = ""
	for num in range(len(probsDict)):
		if probsDict[num][0] > maxprob:
			maxprob = probsDict[num][0]
			maxkmer = probsDict[num][1]
	return maxkmer

"""Example:

tietestprofile =     {"A":[1.0, 1.0, 1.0],"C":[0.0, 0.0, 0.0],"G":[0.0, 0.0, 0.0],"T":[0.0, 0.0, 0.0]} #TEST SUCCESSFUL
tietesttext = "AACCGGTT"  #3-mer 

ProfileMostProbablePattern(tietesttext, 3, tietestprofile) = "AAC"

profile = {"A" : [0.2, 0.2, 0.3, 0.2, 0.3], "C" : [0.4, 0.3, 0.1, 0.5, 0.1], "G" : [0.3, 0.3, 0.5, 0.2, 0.4], "T" : [0.1, 0.2, 0.1, 0.1, 0.2]}
Text = "TTACCATGGGACCGCTGACTGATTTCTGGCGTCAGCGTGATGCTGGTGTGGATGACATTCCGGTGCGCTTTGTAAGCAGAGTTTA" #5-mer

ProfileMostProbablePattern(Text, 5, profile) = "CAGCG"
""" 
######################################################################### GreedyMotifSearch ################################################################################################# 

#The greedy algorithm to find the most conserved motif (i.e motif with the least score)
#t is the no of strings in Dna
#Starting with the first k-mer of the first string, the algorithm finds the profile-most-probable k-mer in the next strings and hence calculates the score.
#However the algorithm does not work as if the k-mer is not present in the second, the first string, as all have prob 0, is chosen.
def GreedyMotifSearch(Dna, k, t):
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

"""Example: 

DNA = ["GACCTACGGTTACAACGCAGCAACCGAAGAATATTGGCAA",
       "TCATTATCGATAACGATTCGCCGGAGGCCATTGCCGCACA",
       "GGAGTCTGGTGAAGTGTGGGTTATGGGGCAGACTGGGAAA",
       "GAATCCGATAACTGACACCTGCTCTGGCACCGCTCTCATC",
       "AAGCGCGTAGGCGCGGCTTGGCATCTCGGTGTGTGGCCAA",
       "AATTGAAAGGCGCATCTTACTCTTTTCGCTTAAAATCAAA",
       "GGTATAGCCAGAAAGCGTAGTTAATTTCGGCTCCTGCCAA",
       "TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT"]
GreedyMotifSearch(DNA, 4, 8) = ['CGCA', 'CGCA', 'GGAG', 'GGCA', 'GGCA', 'CGCA', 'GGTA', 'GGCA']
"""


#Pseudocounts eliminate the problem of probability 0 (By Laplace's Rule of Succession)
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



def GreedyMotifSearchWithPseudocounts(Dna, k, t):
	BestMotifs = []
	for strindex in range(0, t):
		BestMotifs.append(Dna[strindex][0:k])
	strlen = len(Dna[0])
	for i in range(strlen-k+1):
		Motifs = []
		Motifs.append(Dna[0][i:i+k])
		for strindex in range(1, t):
			P = ProfileWithPseudocounts(Motifs[0:strindex])
			Motifs.append(ProfileMostProbablePattern(Dna[strindex], k, P))
		if Score(Motifs) < Score(BestMotifs):
			BestMotifs = Motifs
	return BestMotifs

#Dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG"]
#print GreedyMotifSearchWithPseudocounts(Dna, 3, 5)

#Given an arbitrary profile, find the profile-most probable k-mer in each string of DNA
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

"""Example:
Dna = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT"]
Profile = {"A" : [0.8, 0.0, 0.0, 0.2],
		   "C" : [0.0, 0.6, 0.2, 0.0],
		   "G" : [0.2, 0.2, 0.8, 0.0],
		   "T" : [0.0, 0.2, 0.0, 0.8]}
Motifs(Profile, Dna) = ['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT']
"""

def RandomMotifs(Dna, k, t):
	M = len(Dna[0])-k
	RandMots = []
	for strindex in range(t):
		randpos = random.randint(1, (M+1))
		i = randpos - 1
		RandMots.append(Dna[strindex][i:i+k])
	return RandMots

"""Example
RandomMotifs(["ATCC", "CTGG", "GGGA", "GTAC"], 2, 4) = ['CC', 'CT', 'GG', 'AC']
"""

############################################################################### Randomized Algorithm for Motif Finding #####################################################################################

#Randomized Algorithm for Motif Finding: a collection of randomly chosen k-mers from each string is chosen as the starting BestMotifs. By creating a Profile and finding the 
"""
profile-most probable k-mer from each string(Motifs(Profile(motifs)etc.), we get a BestMotifs with lower and lower score, up to a point where it has found the local optimum. 
If RandomisedMotifSearch is run many times, we can choose the best one from many returned BestMotifs.
In particular, only if one of the implanted motifs is chosen as a random motif, then that may eventually lead the algorithm closer, or to the actual set of implanted motifs.
If Dna was random then the profile would be completely uniform. 
But Dna is not random because of the higher probability of the motifs implanted and the profile points toward the implanted IF one of the implanted motifs is chosen by chance. 
That is its weakness as it is guaranteed to reach the best scoring collection of Motifs that can be reached from its starting point (local optimum) but if its randomly chosen starting
point is wrong then it cannot reach the best solution from there, which is why it is run many times.

""" 
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
		if Dict[i][1] < minscore:
			minscore = Dict[i][1]  
			minbestmots = Dict[i][0]
	return minbestmots

"""Example:
k=8 
t=5
Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
 "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
 "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
 "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
 "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

RepeatedRMotifSearch(Dna, k,t, 100) = ['GGGTGTTC', 'ATGTGTAA', 'AAGTATAC', 'AGGTGCAC', 'ACGTGCAA']"""


#########################################################################   Gibbs Sampling algorithm   #############################################################################################

#Ouput: Sum of numbers in each column = 1
def Normalize(Probabilities):
	NormDict = {}
	ProbSum = 0
	for key in Probabilities:
		ProbSum += Probabilities[key]
	for key in Probabilities:
		NormDict[key] = Probabilities[key]/ProbSum
	return NormDict

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers occurring in DNA
# Output: chosen k-mer with respect to the values in Probabilities
#uniform(a,b) returns a random floating point number
def WeightedDie(Probabilities):
	randNum = random.uniform(0,1)
	kmer = ""
	Probs = Normalize(Probabilities)
	SortedDict = {}
	RangeDict = {}
	i = 0
	for key in Probs:
		SortedDict[i] = [key, Probs[key]]
		i += 1
	LowRange = 0
	for n in range(0, (len(SortedDict))):
		SDentry = SortedDict[n]
		UpRange = LowRange + SortedDict[n][1]
		RangeDict[n] = [ SDentry[0], [LowRange, UpRange] ]
		LowRange += SDentry[1]
	for n in range(len(RangeDict)):
		if randNum >= RangeDict[n][1][0] and randNum < RangeDict[n][1][1]:
			kmer = RangeDict[n][0]
			break
	return kmer

def ProfileGeneratedString(Text, profile, k):
	lTxt = len(Text)
	probabilities = {}
	for i in range(0,lTxt-k+1):
		probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
	probabilities = Normalize(probabilities)
	return WeightedDie(probabilities)

"""Example:
profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
ProfileGeneratedString("AAACCCAAACCC", profile, 2) = CC
"""


# t is the no. of strings in Dna
#Motifs is defined at the start as a collection of random motifs, one from each string. Then one random string is chosen. 
#Motifswithout is basically Motifs without the random motif from the string chosen in the previous step.
#Then a profile with pseudocounts for Motifswithout and an empty dictionary are created. Then the dictionary is filled with each k-mer in the omitted string and its probability 
#according to the profile of the rest of the strings.
#Then the original Motifs entry for the omitted string is replaced by choosing a k-mer using the previously made dictionary and WeightedDie.
#If the new Motifs is better than the old BEstMotifs, then BestMotifs is replaced. 
#This algorithm is better than RandomisedMotifSearch because it doesn't just throw away the entire Motifs after every round of searching.
#It slowly changes one motif at a tme to try and reach the best BestMotifs of all rather than the local optimum- which means that it is more efficient than RAndomisedMotifSearch. 
def GibbsSampler(Dna, k, t, N):
	Motifs = RandomMotifs(Dna, k, t)
	BestMotifs = Motifs
	for r in range(1, N):
		sIndex = random.randint(0, (t-1))
		Motifswithout = []
		for i in range(len(Motifs)):
			if not (i == sIndex):
				Motifswithout.append(Motifs[i])
		profile = ProfileWithPseudocounts(Motifswithout)
		probabilities = {}
		ostr = Dna[sIndex]#omitted string
		for i in range(0, len(Dna[sIndex])-k+1):
			probabilities[ostr[i:i+k]] = Pr(ostr[i:i+k], profile)
		Motifs[sIndex] = WeightedDie(probabilities)
		if Score(Motifs) < Score(BestMotifs):
			BestMotifs = Motifs
	return BestMotifs

def RepeatedGSamplerSearch(Dna, k, t, N, n):
	Dict = {}
	for r in range(n):
		Result =  GibbsSampler(Dna, k , t, N)
		Dict[r] = [Result, Score(Result)]
	minscore = 100000000
	minbestmots = []
	for i in range(len(Dict)):
		#print Dict[i][1], minscore 
		if Dict[i][1] < minscore:
			minscore = Dict[i][1]  
			minbestmots = Dict[i][0]
	return minbestmots

Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
	   "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
	   "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
	   "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
	   "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

"""GibbsSampler(Dna, 8, 5, 100) = ['AAACGGCC', 'CAAGGTGC', 'AAGTATAC', 'TCAGGTGC', 'CCAGCTCC']"""



