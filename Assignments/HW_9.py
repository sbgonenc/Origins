import pandas as pd
from optparse import OptionParser
import matplotlib.pyplot as plt

def data_parser(text, day):
	'''
	parser for tsv files.  drops non-gene ids
	:param text: tsv files
	:return: pandas dataframe
	'''
	genes = []
	tpm_lst = []
	colum = ["gene_id", f"{day}"]
	df = pd.DataFrame(columns= colum)
	with open(text, "r") as fh:
		for lines in fh:

			if not lines.startswith("E"): continue
			words = lines.strip().split()

			gene_id = words[0]
			tx_ids = list(words[1].split(","))
			tpm = words[5]

			genes.append(gene_id)
			tpm_lst.append(tpm)

		df["gene_id"] = genes
		df[f"{day}"] = tpm_lst
		df.set_index("gene_id", inplace=True)

	return df


#if __name__ == "__main__":
#polyA plus RNAseq data for mouse C57BL/6 for liver samples for each day. Each variable contains gene_id and tpm values with days
day11_2 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E11.5/ENCFF523MEO.tsv', 11.5)
day11_1 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E11.5/ENCFF954EHG.tsv', 11.5)
day12_1 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E12.5/ENCFF669WDQ.tsv', 12.5)
day12_2 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E12.5/ENCFF997HBI.tsv', 12.5)
day13_2 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E13.5/ENCFF336VTP.tsv', 13.5)
day13_1 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E13.5/ENCFF615ZTQ.tsv', 13.5)
day14_1 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E14.5/ENCFF432ZGG.tsv', 14.5)
day14_2 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E14.5/ENCFF572OPZ.tsv', 14.5)
day15_2 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E15.5/ENCFF504YJB.tsv', 15.5)
day15_1 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E15.5/ENCFF740SXP.tsv', 15.5)
day16_2 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E16.5/ENCFF512KYX.tsv', 16.5)
day16_1 = data_parser(r'C:\Users\BERK\PycharmProjects\New_Page/Assignments/HW_9/E16.5/ENCFF759PUL.tsv', 16.5)
days_list = [day11_1, day12_1, day13_1, day14_1, day15_1, day16_1]
days_list2 = [day11_2, day12_2, day13_2, day14_2, day15_2, day16_2]


the_DATA = days_list[0].join(days_list[1:], how="outer")
the_DATA2 = days_list2[0].join(days_list2[1:], how="outer")

#the_DATA.to_csv('Assignments/HW_9/merged_data.txt')
#the_DATA2.to_csv('Assignments/HW_9/merged_data2.txt')


if __name__ == "__main__":
	op = OptionParser()
	op.add_option('-g', '--gene_id', help='Enter a valid gene_id(s). Separate them with commas', dest='gene_id')
	op.add_option('-d', '--days', help='Optional. Add a range of day initial:final', dest='days', default='11.5:16.5')


	(options,args) = op.parse_args()


	_days =[11.5,12.5,13.5,14.5, 15.5, 16.5]

	def days_corrector(day=options.days, daylist=None):
		'''
		converts user input from command prompt to codable one
		:param day: command prompt input or default
		:param daylist: indicated _days list above
		:return: daylist for sketcher
		'''
		if daylist is None:
			daylist = _days

		those_were_the_days_my_friend = day.split(':')
		initial_day = float(those_were_the_days_my_friend[0])
		final_day = float(those_were_the_days_my_friend[1])
		index_init = daylist.index(initial_day)
		index_final = daylist.index(final_day)

		return daylist[index_init:index_final+1]

	def sketcher(gene_id=options.gene_id, days=days_corrector(), data=the_DATA):
		'''
		Sketches a line graph
		:param gene_id: user input from command prompt
		:param days: requires days_corrector func, optional user input
		:param data: the parsed and combined pandas dataframe
		:return: fruitless, sketches a line graph
		'''
		gene_ids = list(gene_id.split(','))
		for genes in gene_ids:
			x_values = [day for day in days]
			y_values = list(map(float, [data.loc[genes][f'{x}'] for x in x_values]))
			plt.ylabel(f'TPM Values')
			plt.plot(x_values, y_values, label=f'{genes}', linewidth=2)
			plt.legend()
		plt.show()

	#sketcher(gene_id='ENSMUSG00000000001.4,ENSMUSG00000000031.11,ENSMUSG00000000441.13', days=[11.5,12.5,13.5,14.5, 15.5, 16.5])
	sketcher()