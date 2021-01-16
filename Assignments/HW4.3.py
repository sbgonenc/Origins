genelist = ['GPR139 ', 'YAP1', 'RASGEF1B', 'PAH', 'PLCB2', 'GAPDH', 'SST']
lengthlist = [1613, 2393, 2277, 4122, 4616, 1875, 618]
def gene_dict():
	dictionary = {}
	for idx, genes in enumerate(genelist):
		dictionary[genes] = int(lengthlist[idx])
	return dictionary
#print(gene_dict())
def compare(num):
	mkeys = []
	toplam = 0
	value_input = num
	for keys in gene_dict():
		if gene_dict()[keys] > value_input:
			mkeys.append(keys)
			toplam += gene_dict()[keys]
	return mkeys, toplam
def interface():
	command = input("Please enter a length (Type \'exit\' to quit)\n")
	if command.isnumeric():
		if int(command) <= 0:
			print("Not a valid length!")
			return interface()
		num = int(command)
		genes, sums= compare(num)[0], compare(num)[1]
		if genes == []:
			print("There is no gene which is that long")
			return interface()
		print("Genes whose length is greater than {} are".format(num))
		for each in genes:
			print(each)
		print("The sums of these genes are ", sums)
		return interface()
	elif command.upper() == "EXIT":
		return None
	else:
		print("Couldn't recognize the input.")
		return interface()
interface()