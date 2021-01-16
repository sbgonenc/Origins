seq_input = "ATATGACTAGCCTGCATTACCAATGAAACCGATCCATTTGCATCCAATACGAACACACAAGCTCGCTGAATTGTCCCTCGTTGCACCGAGCTGCATTGAAATCAGGTATGTTTGGCGACATAGGCCCTGAGCTTTTTCGTGTCCAGTATGGTAGGCCCTTCGGTAAGTCCTCGAAGCCATATTTCCTCATTACTCACATCGGAGCTCGGGCACCAGGAGGGTGACGAAGCATTAGAAACACCCCAGTAGCTAACGGGTCGTGTCTCGTATATCGCTCGAATAAGAGAGGGTAGGCTGGACGACACTGGATTACCGTGGCGAGATTCTTGAATTGCAACCACGCGAGCCAGTATGGGGTGAAGTACAAGGCCGGTAAAAGGGCCATCTTTGCCGCGCAAACAGTAGTTGTCCAAGACACGCGCCCGTTCCTCCGAAAGTTATTATCAGCAGCCGTTAGTGCTAAGTCCTCACTGCTAGGGAGCCGTATATCCTAGCTACCCGGAAACCTCTATCCCTGCCTGCCACCGCCTCGATATAGAGTCAGCCTCACGCAGTTTGTATGAGATAGTAACCGACACTTGCTCCGGAGATGCAACGCTGACTGTTGGGGAGAAGCAGATCTATGAGGGTGAAGGTCCCCAGGTACGCGACAAAATCAATGTCCGCTATGAAGTGATAATTGTCAGTCGCGATTTCAGGATGGTTGCGTCACTCCAATTGTGGGCGCGAGCATTGGAGAGAGTTTACAGCCTATCATGTCAGCAGGGTGGAGTCTGCCGGTTCGCTCCTCATCCTATCCACACGAGTCTCGGGTCCTGGGCGGGGGATGGAAATCTTAAAACTATGCTTTACGTAACGTTACAGACTGTTACGTATGGGGCCTGATTATCGACTCACGTATGCAGGCGCTGGAGA".upper()
#seq_input = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT".upper()
file = open("exampleseqs.txt", "r")
#file.write(">Rosalind_6404\nCCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG\n>Rosalind_5959\nCCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC\n>Rosalind_0808\nCCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")
file.close()

def contenter(seq):
	c = 0
	l = 0
	for letter in seq:
		l += 1
		if letter == "G" or letter == "C":
			c += 1
	rv = c/l
	return rv*100
#print(contenter(seq_input))

def parser():
	seq_id = []
	seq_list = []
	seqss = []
	seq_str = ""
	seq_con = {}
	with open("rosalind_gc.txt") as fh:
		for line in fh:
			if line.startswith(">"):
				seq_id.append(line.strip())
				seq_str = str("".join(seq_list))
				seqss.append(seq_str)
				seq_list.clear()
				continue
			else:
				seq_list.append(line.strip())
		seq_str = str("".join(seq_list))
		seqss.append(seq_str)
		for p, v in enumerate(seq_id):
			seq_con[v] = seqss[p+1]
	return seq_con
print(parser())
def gccontent():
	seq_dict = parser()
	gc_list = []
	gc_dic = {}
	c = 0
	for a in seq_dict.keys():
		c = 0
		for letter in seq_dict[a]:
			if letter == "G" or letter == "C":
				c += 1
		gc_list.append(c*100/len(seq_dict[a]))
		gc_dic[a] = c*100/len(seq_dict[a])
	return gc_dic
print(gccontent())

def maxgetter():
	rdict = gccontent()
	max_value = max(rdict.values())
	for r in rdict.keys():
		if max_value == rdict[r]:
			return r, rdict[r]
print(maxgetter())