def parser():
	line_lst = []
	with open("rosalind_gc.txt", "r") as filehandler:
		for line in filehandler:
			if line.startswith(">"): continue
			line_lst.append(line.strip())
		rv = "".join(line_lst)
	return rv
seq_input = parser()

print(seq_input)


def comp_seq(DNA):
	comp_seq = []
	for letter in DNA:
		if letter == "A":
			comp_seq.append("T")
		elif letter == "T":
			comp_seq.append("A")
		elif letter == "G":
			comp_seq.append("C")
		elif letter == "C":
			comp_seq.append("G")
	comp_seq.reverse()
	c_seq = "".join(comp_seq)
	return c_seq

def reader():
	comseq_str = "".join(comp_seq(seq_input))
	for position, letter in enumerate(seq_input):
			for i in range(9):
				seq_str = seq_input[position:position + 4 + i]
				c_str = comseq_str[len(comseq_str) - position- len(seq_str): len(comseq_str) - position- len(seq_str) + 4 + i]
				if seq_str == c_str:
					print(position+1, len(seq_str))
#reader()
