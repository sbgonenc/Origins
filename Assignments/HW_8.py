'''
Using pandas library do the following:

Create a pandas dataframe using the information in "ensembl_kisa.txt" file. 1
Change the index to the transcript ID. 2
Add three new columns for cds_seq, aa_seq and tx_length 3
Reformat the exon_starts and exon_fields by changing the commas to tabs 4
Fill tx_length field for all transcripts  5
For each transcript, if it is a coding transcript, fill in the cds_seq and aa_seq fields using the class you wrote in
your previous assignment. 6

After you completed the above tasks, print the IDs and gene names of longest and shortest transcripts on the screen 7

'''
import pandas as pd
import time
from Assignments.HW_8 import Transcript_Class
import concurrent.futures

t1 = time.perf_counter()

pd.set_option('display.max_columns', 18)

tc = Transcript_Class.Transcript

tx_file = "Assignments/ensembl_kisa.txt"

tx_length_lst = [x for x in Transcript_Class.output("tx_length", tx_file, filter= False)] #check 5
cds_length_lst = [x for x in Transcript_Class.output("cds_length", tx_file, filter= False)]

dataframe = pd.read_table("Assignments/ensembl_kisa.txt", sep = "\t", header= 0) # check 1
dataframe.rename(columns= {"name": "tx_id"}, inplace=True)
dataframe.rename(columns= {"name2": "gene_id"}, inplace=True)
dataframe.set_index("tx_id", inplace= True)  # check 2

dataframe["tx_length"] = tx_length_lst  # check3 / check5
dataframe["cds_length"] = cds_length_lst   # I wanted to include
dataframe["cds_seq"] = None 	# check 3
dataframe["aa_seq"] = None		# check 3
dataframe["exonStarts"] = dataframe["exonStarts"].apply(lambda x: x.replace(",", "\t")) #check 4
dataframe["exonEnds"] = dataframe["exonEnds"].apply(lambda x: x.replace(",", "\t"))		#check4

#dataframe.drop(columns=["#bin", "cdsStartStat", "cdsEndStat", "exonFrames", "score", "txStart", "txEnd", "cdsStart",	"cdsEnd",	"exonCount",	"exonStarts",	"exonEnds"], inplace= True)
	# did not want to look my dataframe messy



def again_a_parser(text):
	'''
	parses the sequence including txt file
	:param text: sequences text
	:return: dictionary {tx_id: [cds_seq, aa_seq]}
	'''

	seq_dictionary = {}
	with open(text, "r") as fh:

		for line in fh:
			if line.startswith("Tx_id"): continue
			words = line.strip().split()
			tx_id = words[0]
			cds_seq = words[2]
			pro_seq = words[3]
			seq_dictionary[tx_id] = [cds_seq, pro_seq]

	return seq_dictionary

sqdict = again_a_parser("Assignments/HW_8/ensembl19_kisa_seq.txt")




def seq_iterator(seq_dictionary, dataframe):  #Check 6
	'''
	to fill in sequences to corresponding tx_id
	:param seq_dictionary: output from "again_a_parser" function
	:param dataframe:  the dataframe
	:return: filled in dataframe
	'''
	for tx_id in seq_dictionary:

		dataframe.loc[(tx_id, "cds_seq")] = seq_dictionary[tx_id][0]
		dataframe.loc[(tx_id, "aa_seq")] = seq_dictionary[tx_id][1]

	return dataframe



df_appended = seq_iterator(sqdict, dataframe)

df_appended.to_csv("Assignments/HW_8/kisa_filled(all_data).txt", sep= ";", index= True, header= True) # I wanted to record the appended dataframe


def mandm(df):  #check 7
	'''
	maximum and minimum transcript length finder
	## I give cds_length output also, since tx_length changes dramatically from one database to other##
	:param df: appended dataframe
	:return: print(tx_id, gene_id, tx_length, cds_length)
	'''

	tx_lengths = list(df["tx_length"])
	mx_tx_length = max(tx_lengths)
	mn_tx_length = min(tx_lengths)

	max_filter = (df["tx_length"] == mx_tx_length)
	print(df.loc[max_filter, ["gene_id", "tx_length", "cds_length"]])

	min_filter = (df["tx_length"] == mn_tx_length)
	print(df.loc[min_filter, ["gene_id", "tx_length", "cds_length"]])

	#while mn_tx_length == 0:
	#	tx_lengths.remove(0)
	#	mn_tx_length = min(tx_lengths)

	#if mn_tx_length != 0:

	#	min_filter = (df["tx_length"] == mn_tx_length)
	#	print(df.loc[min_filter, ["gene_id", "tx_length", "cds_length"]])



mandm(df_appended)



t2 = time.perf_counter()
print(f'Finished in {round(t2-t1)} seconds')

