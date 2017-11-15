#/usr/bin/env python2.7
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
			if nt.isalpha() :
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
			elif nt.isalpha()!=True:
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
		if i.isalpha() or i in [".", "|", "$", "^"]:
	 		self.reverse_str += i
	 		regex_subseq =regex_subseq[:-1]
 			regex_seq_finder.regex_reverse(self, regex_subseq)
		if i == '}' :
			# todoreverse group1
			regex_group = re.compile(r".*(\(.*?\))(\{.*?\})$")
			match = regex_group.match(regex_subseq)
			if match:
				group2 = match.group(2)
				group1 = match.group(1)
			regex_group_letters = re.compile(r".*(\.|[A-Z]{1})(\{.*?\})$")
			match2 = regex_group_letters.match(regex_subseq)
			if match2:
				group2 = match2.group(2)
				group1 = match2.group(1)
			regex_group_letters = re.compile(r".*(\[.*\])(\{.*?\})$")
			match3 = regex_group_letters.match(regex_subseq)
			if match3:
				group2 = match3.group(2)
				group1 = match3.group(1)
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
			regex_group_bracket = re.compile(r".*(\[.*?\])$")
			match_bracket = regex_group_bracket.match(regex_subseq)
			if match_bracket:
				group = match_bracket.group(1)
				self.reverse_str += group
				tmp_regex_grp = re.escape(group)+'$'
				regex_subseq = re.sub(tmp_regex_grp,'',regex_subseq)
				regex_seq_finder.regex_reverse(self, regex_subseq)

	
		if i == '*' :
			regex_group_star = re.compile(r".*(\(.*?\))(\*|\+)$")
			match_star1 = regex_group_star.match(regex_subseq)
			if match_star1:
				group2 = match_star1.group(2)
				group1 = match_star1.group(1)
			regex_group_star2 = re.compile(r".*(\.|[A-Z]{1})(\*)$")
			match_star2 = regex_group_star2.match(regex_subseq)
			if match_star2:
				group2 = match_star2.group(2)
				group1 = match_star2.group(1)
			regex_group_star3 = re.compile(r".*(\[.*?\])(\*)$")
			match_star3 = regex_group_star3.match(regex_subseq)
			if match_star3:
				group2 = match_star3.group(2)
				group1 = match_star3.group(1)

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

	def create_pattern(self,regexp, nuctype, IUPAC):
		"""
		a function to transform iupac value and compress regexp to a
		re object
		:param regexp: regular expression to compress
		:param nuctype: type of nucleic acid
		:param IUPAC: use or not of IUPAC code
		:return: return compressed regular expression
		"""
		pattern = self.use_iupac(regexp, nuctype)
		compiled_pattern =re.compile(pattern)
		return compiled_pattern


	def use_iupac(self, regexp, nuctype):
		"""

		:param regexp: the regular expression to modify IUPAC
		:param nuctype: type of nucleic acid
		:return: regular expression with iupac data modified
		"""
		pattern=""
		for nt in regexp:

			# ipdb.set_trace()
			if nt.isalpha():
				if nuctype == "DNA":
					value = IUPACData.ambiguous_dna_values[nt]
					if len(value) == 1:
						pattern += value
					else:
						pattern += '[%s]' % value
				if nuctype == 'RNA':
					value = IUPACData.ambiguous_rna_values[nt]
					if len(value) == 1:
						pattern += value
					else:
						pattern += '[%s]' % value
			elif nt.isalpha() != True:
				pattern += nt
		return pattern

	def find_subseq(self, sequence, regex,compiled, IUPAC, number_of_match, position_of_match, match, nuctype="DNA", overlap=False):
		"""
		a function to find a subsequence in a sequence based on a regex. This function
		:param sequence:the sequence where to search the regular expression
		:param regex: the regular expression to find
		:param compiled : boolean value to say if the regex is compiled or not
		:param IUPAC : usage of iupac code
		:param number_of_match: a param to set the method to return the number of match between regular expression and sequence
		:param position _of_match: a param to set the method for returning the matching position
		:param match: a param to set the method for returning if there is a match or not
		:param overlap: a param for choosing if it can have an overlap between 2 match
		:type sequence: string
		:type regex: string
		:type IUPAC : boolean
		:type number_of_match: boolean
		:type position_of_match: boolean
		:type match: boolean
		:type overlap: boolean

		"""

		#ipdb.set_trace()
		self.sequence = sequence
		self.nuctype = nuctype

		pattern = ""
		cpt=0
		if compiled !=True:
			self.regex_subseq=regex
			pattern = self.create_pattern(regex,self.nuctype,IUPAC)

		else:
			pattern=regex
			self.regex_subseq = regex.pattern

		#ipdb.set_trace()
		self.find_subseq_result.append(pattern.pattern)
		if number_of_match == True:
			self.find_subseq_result.append(len(pattern.findall(self.sequence,overlap)))
			return self.find_subseq_result

		if match == True :
			match_result = False
			#ipdb.set_trace()
			if pattern.search(self.sequence) :
				match_result =True
			self.find_subseq_result.append(match_result)
			return self.find_subseq_result

		if position_of_match == True:
			matches = pattern.finditer(self.sequence)
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

