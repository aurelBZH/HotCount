#/usr/bin/python
import re
from Bio.Data import IUPACData
import ipdb

class regex_seq_finder(object):
	"""docstring for regex_seq_finder"""
	def __init__(self, arg):
		self.arg = arg
		self.sequence = ""
		self.regex_subseq = ""
		self.nuctype = "DNA"
		self.find_subseq_result = []
		

	
	def find_reverse_complement(self, regex_subseq, type):
		self.regex_subseq = regex_subseq


	def find_subseq(self, sequence, regex, number_of_match, position_of_match, nuctype="DNA", overlap=False):
		self.sequence = sequence	
		self.nuctype = nuctype
		self.regex_subseq = regex	
		pattern = ""
		cpt=0 
		for nt in regex:
			cpt = cpt+1
			# ipdb.set_trace()
			if re.match(r"[A-Z]", nt) :
				if nuctype == "DNA": 
					value = IUPACData.ambiguous_dna_values[nt] 
					if len(value) == 1: 
						pattern += value
					else: 
						pattern += '[%s]' %value 
			elif re.match(r"[^A-Z]", nt):
				pattern += nt 
		self.find_subseq_result.append(pattern)
		if number_of_match ==True :
			number_subseq = 0 
			if len(re.findall(pattern, self.sequence,overlap))>0 :
				number_subseq+=1 
			self.find_subseq_result.append(number_subseq) 
			return self.find_subseq_result 

		if position_of_match == True:
			matches=re.finditer(pattern, self.sequence)
			for match in matches :
				self.find_subseq_result.append(match.start())
			return self.find_subseq_result


if __name__ == '__main__':
	print(regex_seq_finder(object).find_subseq("ATCTTTTTATTTCGCGCGGGGAAA",r"AW{1,10}(CG){1,10}", True, False ))