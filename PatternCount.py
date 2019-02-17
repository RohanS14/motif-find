#Look for the origin of replication and DNAa boxes


#Lists the most frequent patterns of length k in a given Text

def remove_duplicates(List):
    SingleList = []
    for i in List:

        if i not in SingleList:
            SingleList.append(i)
    print SingleList   

def FrequentPatterns(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsndup = remove_duplicates(FrequentPatterns)
    return FrequentPatternsndup 


def CountDict(Text, k):
    Count = {} 
    for i in range(len(Text) - k + 1):
        Pattern = Text[i : i + k]
        Count[i] = PatternCount(Pattern, Text)
    return Count


def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

#words = FrequentPatterns()
#finds the starting positions of a Pattern in Text

def PatternMatching(Pattern,Genome):
	positions = []
	for i in range(len(Genome)-len(Pattern)+1):
		if Genome[i:i+len(Pattern)] == Pattern:
			positions.append(i)
	print positions

#Creates a dictionary with a key for each position on the Genome to keep track of C and G occurences- reduced by 1 if C and increased by 1 if G
#Plotting a graph of the result of this functions(Skew Diagram) shows lines with a positive gradient - which represent occurences of G and lines with a negative gradient - C.
#Due to deamination this can help locate the lagging and leading strands.
def Skew(Genome):
	skew = {}
	skew[0] = 0
	for i in range(len(Genome)):
		if Genome[i] == "C":
			skew[i+1] = skew[i]-1
		elif Genome[i] == "G":
			skew[i+1] = skew[i]+1
		else:
			skew[i+1] = skew[i]
	return skew

#Skew("CATGGGCATCGGCCATACGCC")
def MinSkew(Genome):
	positions = []
	skew = Skew(Genome)
	skewmin = max(skew.values())
	for i in range(len(skew)):
		#skew dictionary iterable as keys are consecutive integers
		if skew[i] == skewmin:
			positions.append(i)
	return positions

#print MinSkew("CATTCCAGTACTTCGATGATGGCGTGAAGA")
#print Skew("CATTCCAGTACTTCGATGATGGCGTGAAGA")



#Finds the reverse complement of a Pattern

def reverse(text):
    newstring = ""
    for i in text:
        req_letter = i
        newstring =  req_letter + newstring 
    return newstring


def reverseComplement(text):
	DNA_Strand = text.upper()
	RevCom = ""
	for i in DNA_Strand:
		if i == "A":
			RevCom = RevCom + "T"
		elif i == "T":
			RevCom = RevCom + "A"
		elif i == "C":
			RevCom = RevCom + "G"
		elif i == "G":
			RevCom = RevCom + "C"

	return reverse(RevCom)

"""print reverseComplement( "GCTAGCT")
	
ReverseComplement = reverseComplement("")		
print ReverseComplement"""


#These were used to finally find DNA-A boxes
#Finds the number of differences between two strings of equal length
def HammingDistance(p,q):
	str1 = p
	str2 = q
	hdist = 0
	for n in range(len(str1)):
		if str1[n] != str2[n]:
			hdist += 1
	return hdist

#print HammingDistance("TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC","GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA") 

def ApproximatePatternMatching(Pattern, Text, d):
	positions = []
	for i in range(len(Text)-len(Pattern)+1):
		if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
			positions.append(i)
	return positions

#print ApproximatePatternMatching("ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3)

def ApproximatePatternCount(Pattern, Text, d):
	count = 0
	for i in range(len(Text)-len(Pattern)+1):
		if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
			count = count + 1
	return count

print ApproximatePatternCount("GTGCCG", "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", 3)
