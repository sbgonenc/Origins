sequence = "CGGCCACGGTAGCCTAGGGGTGCCAGGTGAACCCCCATGGCGAAGAGCTGCGGACGGGTAGAAGCCTGGCATTACGATATCGATGTATTGTTAACGCATCGGGACCGATCGAGATACAATGGAATACAGGTTCCTCGGGCCATTAATGCTTCCCCCGTGGTTTATAACACTCCTAAACTACTTCTGGATGTCAGTGTCTCGTGATTGTTATACACACATGTACCGTGAGCATATGAACTACGGCTTGATCGACGGGCACTTTCCGTCGCGGAATAGCAGTCACGTCTGATAAAGAGAGATCCGTCGATGTTCAAGGGTTCTAATCACTATGGTGCTTGGCGCCCAGGATTACGCTCTCAACGTATCGACGTCCGACATTATGGGCAGCCCCATGACTGGGAAGAAGACATAACACGAAATTGGGAATTATGGGGCTGCTGTCATGCAAGCGGTAGGCCTCGGATTAGTCGCAGATGAGACTTCTTGAGTCGATGCTAATGAAGTAATACCACCTTTCAACACAGCCTGAAGTCCACATGAGTCAGGATGGATGACACGTAGGTCAAGGTGACATTACACCGTAGAAGCCAGTTCACATCGTCGACATGATTTTCACAGGCGGCGCCAGCTGCTGTGCACGAGCTAACCGCGGTGGCGTGTACTTACTTCTCGGGGACGCGGGCCTTTCCAAACCGACCCTAGCAAAAGAACACAGAACCTAAGCCTCCAAGCCATTCGGATTCTCAACGATTGCGAGGTGAACACTACAGCATTCCTCGGTCTCGCAGATAGTATTGGATAAGCACCACCGATAGGCGTCTGTCAATGTCCGTGTGCGTCAATTTCAAGCGTAGTCAAGATCGGTGGTGAACGGAAACCATGAGCCCGCTTTGTCGGCCAGACGCGGACAATCAAAGATCCCCCGAGGTGGCGACATTCTCAACGCGGCGCTGGGTCTGCCGAACGTGAGGTCATTTAGAGGATGGCGCCTACC"

def Transcription():
	seq_list = []
	for letter in sequence:
		if letter == "T":
			seq_list.append("U")
		else:
			seq_list.append(letter)
	seq = "".join(seq_list)
	print(seq)

Transcription()