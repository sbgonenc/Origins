from pprint import pprint
import requests
from functools import lru_cache


def fetch_sequence_fromtogows(chrom, start, end, genomeversion='hg19'):
	url = f"http://togows.org/api/ucsc/{genomeversion}/chr{chrom.lower().replace('chr', '')}:{start}-{end}"
	return requests.get(url).content.decode()


def comp_seq(seq_input):
	comp_seq = ""

	for letter in seq_input.upper():

		if letter == "A":
			comp_seq += "T"
		elif letter == "T":
			comp_seq += "A"
		elif letter == "G":
			comp_seq += "C"
		elif letter == "C":
			comp_seq += "G"
		elif letter == "N":
			comp_seq += "N"

		#comp_seq.reverse()
		#c_seq = "".join(comp_seq)

	return comp_seq


class Transcript:

	num_of_tx = 0

	def __init__(self, gene_id, tx_id, chrom, strand, exon_C, exon_S, exon_E, tx_S, tx_E, cds_S, cds_E):
		self.gene_id = gene_id
		self.tx_id = tx_id
		self.chrom = chrom
		self.strand = strand
		self.exon_C = exon_C
		self.exon_S = exon_S
		self.exon_E = exon_E
		self.tx_S = tx_S
		self.tx_E = tx_E
		self.cds_S = cds_S
		self.cds_E = cds_E

		self.cds_seq = ""
		self.cds_coors = []

		Transcript.num_of_tx += 1


	@lru_cache(maxsize= 1000)
	def tx_length(self):
		#Takes exon starts and ends, returns sum(end-start)
		length = 0

		for i in range(self.exon_C):
			length += self.exon_E[i] - self.exon_S[i] # +1

		return length

	def cds_length2(self):
		'''
		More reliable way to calculate cds length. Also positions cds_start-exon_end-exon_start-cds_end
		:return: cds length
		'''
		cds_start = self.cds_S
		cds_end = self.cds_E
		length = 0

		aralik = []          # to obtain in which exon cds_start, cds_end reside
		cds_cont_coors = []  # for obtaining cds_start, exon_end, exon_start, cds_end coordinates

		for i in range(self.exon_C):
			exon_start = self.exon_S[i]
			exon_end = self.exon_E[i]

			if cds_start > exon_end: continue

			elif exon_start <= cds_start <= exon_end:
				aralik.append(i)

				if exon_start <= cds_end <= exon_end:

					cds_cont_coors.append([cds_start, cds_end])
					length += cds_end - cds_start
					aralik.append(i)
					break

				else:
					cds_cont_coors.append([cds_start, exon_end])
					length += exon_end - cds_start

			elif exon_start <= cds_end <= exon_end:

				aralik.append(i)
				length += cds_end - exon_start

				cds_cont_coors.append([exon_start, cds_end])
				break

			else:

				length += exon_end - exon_start

				cds_cont_coors.append([exon_start, exon_end])


		#length += 1 --> note that +1 excluded, since the coordinate difference gives exact length (checked from ensembl)

		self.cds_coors = cds_cont_coors

		return length #, aralik, cds_cont_coors

	def cds_sequence(self):

		coordinates = self.cds_coors
		_cds_seq = ""

		for couple in coordinates:

			_cds_seq += fetch_sequence_fromtogows(self.chrom, couple[0]+1, couple[1])

		if self.strand == "+":

			self.cds_seq = _cds_seq

			return _cds_seq

		elif self.strand == "-":

			_cds_seq = comp_seq(_cds_seq)
			self.cds_seq = _cds_seq[::-1]

			return _cds_seq[::-1]

	def protein_seq(self):

		seq = self.cds_seq

		rv = ''
		table = {
			'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
			'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
			'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
			'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
			'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
			'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
			'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
			'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
			'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
			'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
			'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
			'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
			'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
			'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
			'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
			'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
		}

		for idx in range(int(len(seq) / 3)):
			codon = seq[idx * 3:idx * 3 + 3]

			if table[codon] == '_':
				break
			rv += table[codon]

		return rv


def parse(text):
	#To parse ensenmbl data
	statement_42 = [42*x for x in range(13)]
	with open(text) as fh:
		dic={}

		for paragraph, line in enumerate(fh):

			if line.startswith("#"): continue

			boluk=line.strip().split()
			tx_id=boluk[1]
			ts_start = int(boluk[4])
			ts_end = int(boluk[5])
			cds_start = int(boluk[6])
			cds_end = int(boluk[7])
			exon_count = int(boluk[8])
			exon_start = list(map(int, boluk[9].split(",")[:-1]))
			exon_end = list(map(int, boluk[10].split(",")[:-1]))
			strand = boluk[3]
			chrom = boluk[2]
			gene_id = boluk[12]

			if cds_start == cds_end: continue #Filter for 0 transcripts
			if cds_start == ts_start: continue #Filter for pseudogenes/nonsenses
			if cds_end == ts_end: continue #Filter for pseudogenes/nonsenses
			if cds_end <= cds_start + 12: continue #
			#if paragraph not in statement_42: continue

			dic[tx_id] = Transcript(gene_id, tx_id, chrom, strand, exon_count, exon_start, exon_end, ts_start,ts_end, cds_start, cds_end)

	return dic

ensembl = parse("Assignments/ensembl_hg19.txt")


def seq_writer(tx_seqs):
#Gathers sequence from togows, using Transcript class

	with open(tx_seqs, "w") as fh:
		fh.write("Tx_id\tchr\tstrand\ttx_length\tcds_length\tcoding_sequence\ttranslation\n")

		for e in ensembl:
			fh.write(f"{ensembl[e].tx_id}\t{ensembl[e].chrom}\t{ensembl[e].strand}\t{ensembl[e].tx_length()}\t{ensembl[e].cds_length2()}\t{ensembl[e].cds_sequence()}\t{ensembl[e].protein_seq()}\n")

#seq_writer(tx_seqs = "deneme_tahtasi.txt")

cds_length_sum = 0
tx_length_sum = 0
for e in ensembl:
	cds_length_sum += ensembl[e].cds_length2()
	tx_length_sum += ensembl[e].tx_length()

avg_cds = cds_length_sum/Transcript.num_of_tx
avg_tx = tx_length_sum/Transcript.num_of_tx
print("Avg cds length:", avg_cds, "\nAvg tx length:", avg_tx, "\nratios:", avg_cds/avg_tx, Transcript.num_of_tx)