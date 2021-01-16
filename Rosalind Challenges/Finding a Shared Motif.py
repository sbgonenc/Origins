def fasta_parser(txt):
	seq_id = []
	seq_list = []
	seqss = []
	seq_str = ""
	seq_con = {}
	with open(txt) as fh:
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

print(fasta_parser("Assignments/seqs"))

s = """>Rosalind_1
GATTACA
>Rosalind_2
TAGACCA
>Rosalind_3
ATACA""".split(">")[1:]
for i in range(len(s)):
    s[i] = s[i].replace("\n", '')
    while s[i][0] not in "ACGT":
        s[i] = s[i][1:]
# ^^^^^^^^^^^^^ all of that to format in FAST in array

#Get shortest of DNA strings
index = s.index(min(s, key=len))

motif = 'A'
shortest = s[index]

#cycle over the DNA string letters
for i in range(len(shortest)):
    n = 0
    present = True
    while present:
            #cycle inside over all other DNA strings and if it's present in there considered a motif and length gets increased by 1
        for each in s:
            if shortest[i:i+n] not in each or n>1000:
                present = False
                break
        if present:
            motif = max(shortest[i:i+n], motif, key=len)
        n += 1
print (motif)