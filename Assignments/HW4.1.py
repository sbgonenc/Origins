file = "ensembl_hg19.txt"

def Avg_5UTR():
	Upstream = []
	Downstream = []
	with open(file) as filehandler:
		for line in filehandler:
			if line.startswith("#"): continue
			words = line.strip().split("\t")
			if words[3] == "+":						# strand <- word[3]
				UTR5 = int(words[6]) - int(words[4]) # cdsStart <- word[6], txStart <- word[4]
				Upstream.append(UTR5)
			elif words[3] == "-":
				UTR5 = int(words[5]) - int(words[7])  # cdsEnd <- word[7], txEnd <- word[5]
				Downstream.append(UTR5)
		#up_avg = sum(Upstream)/len(Upstream)
		#down_avg = sum(Downstream)/len(Downstream)
		Total = sum(Upstream)+ sum(Downstream)
		total_l = len(Upstream) + len(Downstream)
	return Total/total_l
#print(Avg_5UTR())
def Avg_3UTR():
	Upstream = []
	Downstream = []
	with open(file) as filehandler:
		for line in filehandler:
			if line.startswith("#"): continue
			words = line.strip().split("\t")
			if words[3] == "+":						# strand <- word[3]
				UTR3 = int(words[5]) - int(words[7]) # cdsStart <- word[6], txStart <- word[4]
				Upstream.append(UTR3)
			elif words[3] == "-":
				UTR3 = int(words[6]) - int(words[4])  # cdsEnd <- word[7], txEnd <- word[5]
				Downstream.append(UTR3)
		#up_avg = sum(Upstream)/len(Upstream)
		#down_avg = sum(Downstream)/len(Downstream)
		Total = sum(Upstream)+ sum(Downstream)
		total_l = len(Upstream) + len(Downstream)
	return Total/total_l
#print(Avg_3UTR())