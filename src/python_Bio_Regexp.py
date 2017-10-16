#/usr/bin/python
# coding: utf-8 
from __future__ import print_function
import re
from Bio.Data import IUPACData
import ipdb


class regex_seq_finder(object):
	"""
	a class to find and count mutation based on regular expression
	"""
	def __init__(self):
		"""
		init method
		"""
		self.sequence = ""
		self.regex_subseq = ""
		self.nuctype = "DNA"
		self.find_subseq_result = []
		self.complement = ""
		self.reverse_str=""


	def regex_complement(self, regex_subseq, nuctype="DNA"):
		""" 
		a class to find the complement of a regular expression working for DNA and RNA
		:param regex_subseq: the content af the regular expression
		:param nuctype: type of nucleic acid
		:type regex_subseq: string
		:type nuctype: string
		"""
		self.regex_subseq = regex_subseq
		pattern = ""
		for nt in self.regex_subseq:
			if re.match(r"[A-Z]", nt) :	
				if nuctype == "DNA": 
					complement_value = IUPACData.ambiguous_dna_complement[nt]
					value = IUPACData.ambiguous_dna_values[complement_value]
					
				if nuctype == "RNA":
					complement_value = IUPACData.ambiguous_rna_complement[nt]
					value = IUPACData.ambiguous_rna_values[complement_value]
				if len(value) == 1: 
					pattern += value
				else: 
					pattern += '[%s]' %value 
			elif re.match(r"[^A-Z]", nt):
				pattern += nt
		self.complement=pattern	
		return self.complement


	def regex_reverse(self, regex_subseq):
		"""
		a function to make the reverse of a regular expression
		:param regex_subseq: the content of a regular expression
		:type regex_subseq: string
		"""
		if len(regex_subseq) == 0:
			return self.reverse_str
		
		i=regex_subseq[-1]
		if re.match(r"[A-Z]",i) or i in [".", "|", "$", "^"]:
	 		self.reverse_str += i
	 		regex_subseq =regex_subseq[:-1]
 			regex_seq_finder.regex_reverse(self, regex_subseq)
		if i == '}' :
			# todoreverse group1
			if re.match(r".*(\(.*?\))(\{.*?\})$",regex_subseq):
				group2 = re.match(r".*(\(.*?\))(\{.*?\})$",regex_subseq).group(2)
				group1 = re.match(r".*(\(.*?\))(\{.*?\})$",regex_subseq).group(1)
			elif re.match(r".*(\.|[A-Z]{1})(\{.*?\})$",regex_subseq):
				group2 = re.match(r".*(\.|[A-Z]{1})(\{.*?\})$",regex_subseq).group(2)
				group1 = re.match(r".*(\.|[A-Z]{1})(\{.*?\})$",regex_subseq).group(1)

			elif re.match(r".*(\[.*\])(\{.*?\})$",regex_subseq):
				group2 = re.match(r".*(\[.*?\])(\{.*?\})$",regex_subseq).group(2)
				group1 = re.match(r".*(\[.*?\])(\{.*?\})$",regex_subseq).group(1)
			tmp_regex_grp2 = re.escape(group2)+'$'
			tmp_regex_grp1 = re.escape(group1)+'$'
			regex_subseq = re.sub(tmp_regex_grp2,'',regex_subseq)
			regex_subseq = re.sub(tmp_regex_grp1+"$",'',regex_subseq)
		 	self.reverse_str = regex_seq_finder.regex_reverse(self, group1)+group2 
		 	regex_seq_finder.regex_reverse(self, regex_subseq)

	 	if i == ')' :
	 		self.reverse_str += "("
	 		regex_subseq = regex_subseq[:-1]
	 	 	regex_seq_finder.regex_reverse(self, regex_subseq)

		if i == '(':
			self.reverse_str += ")"
			regex_subseq = regex_subseq[:-1]
	 	 	regex_seq_finder.regex_reverse(self, regex_subseq)

		if i == ']' :	
			if re.match(r".*(\[.*?\])$", regex_subseq):
				group = re.match(r".*(\[.*?\])",regex_subseq).group(1)
				self.reverse_str += group
				tmp_regex_grp = re.escape(group)+'$'
				regex_subseq = re.sub(tmp_regex_grp,'',regex_subseq)
				regex_seq_finder.regex_reverse(self, regex_subseq)

	
		if i == '*' :
			if re.match(r".*(\(.*?\))(\*|\+)$",regex_subseq):
				group2 = re.match(r".*(\(.*?\))(\*|\+)$",regex_subseq).group(2)
				group1 = re.match(r".*(\(.*?\))(\*|\+)$",regex_subseq).group(1)

			elif re.match(r".*(\.|[A-Z]{1})(\*)$",regex_subseq):
				group2 = re.match(r".*(\.|[A-Z]{1})(\*|\+)$",regex_subseq).group(2)
				group1 = re.match(r".*(\.|[A-Z]{1})(\*|\+)$",regex_subseq).group(1)

			elif re.match(r".*(\[.*?\])(\*)$",regex_subseq):
				group2 = re.match(r".*(\[.*?\])(\*|\+)$",regex_subseq).group(2)
				group1 = re.match(r".*(\[.*\])(\*|\+)$",regex_subseq).group(1)	

			tmp_regex_grp2 = re.escape(group2)+'$'
			tmp_regex_grp1 = re.escape(group1)+'$'
			regex_subseq = re.sub(tmp_regex_grp2,'',regex_subseq)
			regex_subseq = re.sub(tmp_regex_grp1,'',regex_subseq)
		 	self.reverse_str = regex_seq_finder.regex_reverse(self, group1)+group2 

			regex_seq_finder.regex_reverse(self, regex_subseq)

		return self.reverse_str	

	def verify_regex(self,regex, nuctype = "DNA"):
		"""
		:param regex: the content of a regular expression
		:param nuctype: type of nucleic acid
		:type regex: string
		:type nuctype: string
		"""
		cpt_parenthesis = 0
		cpt_bracket = 0 
		cpt_brace = 0
		if nuctype == "DNA":
			nuc = IUPACData.ambiguous_dna_values
		elif nuctype == "RNA":
			nuc = IUPACData.ambiguous_rna_values	

		for i in regex:
			if i in nuc or  i in [ "{", "}", "(", ")",".", "*", ",", "+", "-", "^","[", "]","$","!", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
				if i  == "(":
					cpt_parenthesis += 1
					continue
				if i ==")":
					cpt_parenthesis -= 1
				if i == "[":
					cpt_bracket +=1
					continue 
				if i == "]":
					cpt_bracket -=1
				if i == "{":
					cpt_brace += 1
					continue
				if i == "}":
					cpt_brace -= 1

				if cpt_brace != 0 :
					if i not in ["0","1","2","3","4","5","6","7","8","9",","]:
						raise Exception("malformed regular expression")
				if cpt_bracket != 0: 
					if i not in nuc and i != "-":
						raise Exception("problem between bracket")
			else : 
				raise Exception ("false symbol in the regular expression")
		if cpt_bracket !=0 or cpt_brace !=0 or cpt_parenthesis !=0:
					raise Exception("anormal number of bracket brace or parenthesis")



	def find_subseq(self, sequence, regex, number_of_match, position_of_match, match, nuctype="DNA", overlap=False):
		"""
		a function to find a subsequence in a sequence based on a regex. This function
		:param sequence:the sequence where to search the regular expression
		:param regex: the regular expression to find
		:param number_of_match: a param to set the method to return the number of match between regular expression and sequence
		:param position _of_match: a param to set the method for returning the matching position
		:param match: a param to set the method for returning if there is a match or not
		:param overlap: a param for choosing if it can have an overlap between 2 match
		:type sequence: string
		:type regex: string
		:type number_of_match: boolean
		:type position_of_match: boolean
		:type match: boolean
		:type overlap: boolean

		"""
		self.sequence = sequence
		self.nuctype = nuctype
		self.regex_subseq = regex	
		pattern = ""
		cpt=0 
		for nt in self.regex_subseq:
		
			# ipdb.set_trace()
			if re.match(r"[A-Z]", nt) :
				if nuctype == "DNA": 
					value = IUPACData.ambiguous_dna_values[nt] 
					if len(value) == 1: 
						pattern += value
					else: 
						pattern += '[%s]' %value 
				if nuctype == 'RNA':
					value = IUPACData.ambiguous_rna_values[nt]
					if len(value) == 1: 
						pattern += value 
					else: 
						pattern += '[%s]' % value
			elif re.match(r"[^A-Z]", nt):
				pattern += nt 
		self.find_subseq_result.append(pattern)

		if number_of_match == True:
			self.find_subseq_result.append(len(re.findall(pattern, self.sequence,overlap)))
			return self.find_subseq_result

		if match == True :
			match_result = False
			if len(re.findall(pattern, self.sequence,overlap))>0 :
				match_result =True
			self.find_subseq_result.append(match_result)
			return self.find_subseq_result 

		if position_of_match == True:
			matches=re.finditer(pattern, self.sequence)
			for match in matches :
				self.find_subseq_result.append(match.start())
			return self.find_subseq_result

	
	def regex_reverse_complement(self, regex, nuctype="DNA"):
		"""
		a simple method to do the reverse and the complement in the same time 
		:param regex: the regular expression to reverse complement
		:param nuctype: the type of nucleic acid
		type regex: string
		type nuctype: string
		"""
		self.verify_regex(regex,nuctype = nuctype)
		reverse_regex = self.regex_reverse(regex)
		
		return self.regex_complement(reverse_regex, nuctype=nuctype)

