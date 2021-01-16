#HW6.1
import numpy as np
from functools import lru_cache

def str_to_int(df):
	#Change string characters to integer in the matrix
	for v, e in enumerate(df):
		for c, f in enumerate(e):
			df[v][c] = float(f)
			df[v][c] = int(df[v][c])
	return df

def regex_conv(m):
	#convert int matrix to a regular expression
	n_matrix = np.array(m)
	tr_matrix = np.transpose(n_matrix)
	consensus = [0 for bah in range(len(tr_matrix))]
	w_list = []
	for i, e in enumerate(tr_matrix):
		for p, a in enumerate(e):
			if a != 0:
				if p == 0: w_list.append("A")
				if p == 1: w_list.append("C")
				if p == 2: w_list.append("G")
				if p == 3: w_list.append("T")
		if len(w_list) > 1:
			w_list.insert(0, "[")
			w_list.append("]")
		consensus[i] = "".join(w_list)
		w_list.clear()
	return "".join(consensus)

@lru_cache(maxsize= 1000)
def filehandler(text = "Assignments/tx_factors.txt"):
	#To parse transcription factors text
	tx_id = []
	empty_list = []
	lst = []
	tx_regex_dic = {}
	with open(text) as fh:
		for lines in fh:
			words = lines.strip().split()
			if lines.startswith(">"):
				tx_id.append(words[1])
				if len(empty_list) == 4:
					int_lst = str_to_int(empty_list)
					lst.append(regex_conv(int_lst))
					empty_list.clear()
					int_lst.clear()
			else: empty_list.append(words)
		if len(empty_list) == 4:
			int_lst = str_to_int(empty_list)
			lst.append(regex_conv(int_lst))
	for i, e in enumerate(tx_id):
		tx_regex_dic[e] = lst[i]
	return tx_regex_dic


#print(filehandler("Assignments/tx_factors.txt"))
dicforpattern = filehandler("Assignments/tx_factors.txt")

#HW6.2

import requests
import pandas as pd

def fetch_sequence_fromtogows(chrom, start, end, genomeversion = 'hg19'):
	url = f"http://togows.org/api/ucsc/{genomeversion}/chr{chrom.lower().replace('chr', '')}:{start}-{end}"
	return requests.get(url).content.decode()

def comp_seq(seq_input):
	comp_seq = []
	for letter in seq_input:
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

@lru_cache(maxsize= 2048)
class transcript():
	def __init__(self, id, ts_start, ts_end, strand, chrom):
		self.id=id
		self.chrom=chrom
		self.strand=strand
		self.ts_start=int(ts_start)
		self.ts_end=int(ts_end)

	def promoter_seq(self):
		if self.strand == "+":
			my_prom_end = self.ts_start
			my_prom_start = my_prom_end - 5001 #5001 olmalı unutma!!!!!!!!!
			promoter = fetch_sequence_fromtogows(self.chrom, my_prom_start, my_prom_end)
		elif self.strand == "-":
			my_prom_end = self.ts_end
			my_prom_start = my_prom_end + 5001 #5001 olmalı unutma!!!!!!!!!!!!!
			promoter = fetch_sequence_fromtogows(self.chrom, my_prom_end, my_prom_start)
			promoter = comp_seq(promoter)
		return promoter

@lru_cache(maxsize= 1000)
def parse(text):
	#To parse ensenmbl dataframe only for chr 22
	with open(text) as sequence:
		dic={}
		for line in sequence:
			if line.startswith("#"): continue
			boluk=line.strip().split()
			id=boluk[1]
			ts_start=boluk[4]
			ts_end=boluk[5]
			#cds_start=boluk[6]
			#cds_end=boluk[7]
			strand=boluk[3]
			chrom=boluk[2]

			if chrom == "chr22":
				dic[id]=transcript(id,ts_start,ts_end, strand, chrom)

	return dic
sdic=parse('Assignments/ensembl_hg19.txt')

#for e in sdic:
#	print(e, sdic[e].chrom, sdic[e].strand, len(sdic[e].promoter_seq()))

#HW6.3
import re
from pprint import pprint

@lru_cache(maxsize = 1024)
def tf_occurance_on_promoter(tf_list=dicforpattern, promoter=sdic):
	'''
	takes two arguments
	:param tf_list: transcription factor patterns
	:param promoter: promoter sequence
	:return: dictionary {gene_id: [{tf: occurance}]}
	'''
	c = 0
	occ_dic = {}
	for sequence in promoter:
		occ_dic.setdefault(sequence, [])

		for tf in tf_list:
			c = 0
			tf_str = tf_list[tf]
			f_iter = re.finditer(tf_str, promoter[sequence].promoter_seq())

			if f_iter != None:
				iter_ls = []
				for x in f_iter:
					iter_ls.append(x)
				c = len(iter_ls)
				#occ_dic[sequence].append({tf:{"Occurance":c, "Spans": iter_ls}})
				occ_dic[sequence].append([tf, c
										  #, iter_ls
											])

		#occ_dic[sequence] = transcript_occurance(sequence, occ_dic[sequence][tf], occ_dic[sequence][tf]["Occurance"], occ_dic[sequence][tf]["Spans"])

	return occ_dic

#toop = tf_occurance_on_promoter()
occ_dic = tf_occurance_on_promoter()
#pprint(occ_dic)
#for e in occ_dic:
#	for b in occ_dic[e]:
#		print(occ_dic[e][0], b)

def co_occurance_finder_each_tx(dictionary):
	'''
	takes a dictionary of tx_id and tf and occurances
	:param dictionary: a dictionary of tx_id and tf and occurances
	:return: co_occurances of tf based on min of two single occurances on each tx
	'''
	co_occurance_dict = {}
	for tx_id in dictionary:
		co_occurance_dict.setdefault(tx_id, [])

		for dicts in dictionary[tx_id]:
			tf1, occurance_number1 = dicts

			for dicts2 in dictionary[tx_id]:
				tf2, occurance_number2 = dicts2

				if tf1 != tf2:
					co_occurance_number = min(occurance_number1, occurance_number2)
					co_occurance_dict[tx_id].append({f"{tf1},{tf2}":co_occurance_number})
				elif tf1 == tf2:
					co_occurance_dict[tx_id].append({f"{tf1},{tf2}": 0})

	return co_occurance_dict

coofet= co_occurance_finder_each_tx(occ_dic)
pprint(coofet)

def sum_of_occurances(co_occurance_dictionary):
	'''
	Takes co_occurance dictionary from tf_occurance_on_promoter, sums the tf occurances
	:param co_occurance_dictionary:
	:return: sum of each tf occurances for each transcripts
	'''

	sum_occurance_dict = {}
	for tx_id in co_occurance_dictionary:

		for tf_list in co_occurance_dictionary[tx_id]:
			#couple = tf_couple.split(",")
			for tf_couple1 in tf_list:
				sum_occurance_dict.setdefault(tf_couple1, 0)

			for tf_couple2, value in tf_list.items():
				sum_occurance_dict[tf_couple2] += value

	return sum_occurance_dict

sums = sum_of_occurances(coofet)
#print("Sums=", sums)


#coofstx = co_occurance_finder_summed_tx()
#pprint(coofstx)



def matrix_converter(co_occurance_dictionary):

	'''
	:param co_occurance_dictionary: Takes a dictionary from  co_occurance_finder
	:return: a co-occurance matrix
	'''

	# Önce dictionary keylerini ayır "," ile.
	rows = []
	#columns = []
	values = []
	for tf_couple, value in co_occurance_dictionary.items():
		couple = tf_couple.split(",")

		if couple[0] == couple[1]:
			values.append(None)
			continue

		if couple[0] not in rows:
			rows.append(couple[0])
		#if couple[1] not in columns: columns.append(couple[1])
		values.append(value)

	data = np.array(values)
	shape = (len(rows), len(rows))
	df = data.reshape(shape)
	dataframe = pd.DataFrame(df,index= rows, columns = rows)

	dataframe.to_csv('tf_co-occurance_chr22.txt', index=True, header=True, sep='\t')
	return dataframe

print(matrix_converter(sums))