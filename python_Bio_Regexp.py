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
		self.complement = ""
		self.reverse_list=""


	def regex_complement(self, regex_subseq, nuctype="DNA"):
		self.regex_subseq = regex_subseq
		pattern = ""
		for nt in self.regex_subseq:
			if re.match(r"[A-Z]", nt) :	
				if nuctype == "DNA": 
					complement_value = IUPACData.ambiguous_dna_complement[nt]
					value = IUPACData.ambiguous_dna_values[complement_value]
					if len(value) == 1: 
						pattern += value
					else: 
						pattern += '[%s]' %value 
			elif re.match(r"[^A-Z]", nt):
				pattern += nt
		self.complement=pattern	
		return self.complement

	def regex_reverse(self, regex_subseq, nuctype="DNA"):
		print regex_subseq+ "muhhh"
		if len(regex_subseq) == 0:
			return self.reverse_list
		regex_list = list(regex_subseq)		
		# for i in reversed(regex_list):
			# ipdb.set_trace()
		i=regex_list[-1]
		print i
		if i in IUPACData.ambiguous_dna_values or i == ".":
			print "test "+i

	 		self.reverse_list += i
				# ipdb.set_trace()

	 		regex_list.pop()
	 		if len(regex_list)>1:
	 			print regex_list
	 			regex_seq_finder.regex_reverse(self, "".join(regex_list), nuctype="DNA")


		if i == '}' :
			# ipdb.set_trace()
			if re.match(r".*(\(.*\))(\{.*\})$","".join(regex_list)):
				# ipdb.set_trace()
				group2 = re.match(r".*(\(.*\))(\{.*\})$","".join(regex_list)).group(2)
				print group2	
				group1 = re.match(r".*(\(.*\))(\{.*\})$","".join(regex_list)).group(1)

			if re.match(r".*([A-Z]{1})(\{.*\})$","".join(regex_list)):
				group2 = re.match(r".*([A-Z]{1})(\{.*\})$","".join(regex_list)).group(2)
				print group2	
				group1 = re.match(r".*([A-Z]{1})(\{.*\})$","".join(regex_list)).group(1)

			if re.match(r".*(\[.*\])(\{.*\})$","".join(regex_list)):
				group2 = re.match(r".*(\[.*\])(\{.*\})$","".join(regex_list)).group(2)
				print group2	
				group1 = re.match(r".*(\[.*\])(\{.*\})$","".join(regex_list)).group(1)	
			for j in group2:
				regex_list.pop()
			for k in group1:
				regex_list.pop()
		 	self.reverse_list += group1+group2 

		 	if len(regex_list)> 1:
		 		# ipdb.set_trace()
		 		# regex_list = self.remove_pattern(r"\{.*\}", regex_list)
			 	regex_seq_finder.regex_reverse(self, "".join(regex_list), nuctype="DNA")

	 	if i == ')' :
	 		group = re.match(r"(\(.*\))$","".join(regex_list)).group(1)
			self.reverse_list += group
			for i in group:
				regex_list.pop()
		if i == ']' :
			if re.match(r".*(\[.* \]$)", "".join(regex_list)):
				group = re.match(r".*(\[.* \])","".join(regex_list)).group(1)
				self.reverse_list += group
			for i in group:
				regex_list.pop()
		# if i == '*' :
		# 	if re.match(r".*(\(.*\))(\*)$","".join(regex_list)):
		# 		# ipdb.set_trace()
		# 		group2 = re.match(r".*(\(.*\))(\*)$","".join(regex_list)).group(2)
		# 		print group2	
		# 		group1 = re.match(r".*(\(.*\))(\*)$","".join(regex_list)).group(1)

		# 	if re.match(r".*([A-Z]{1})(\{.*\})$","".join(regex_list)):
		# 		group2 = re.match(r".*([A-Z]{1})(\*)$","".join(regex_list)).group(2)
		# 		print group2	
		# 		group1 = re.match(r".*([A-Z]{1})(\*)$","".join(regex_list)).group(1)

		# 	if re.match(r".*(\[.*\])(\{.*\})$","".join(regex_list)):
		# 		group2 = re.match(r".*(\[.*\])(\*)$","".join(regex_list)).group(2)
		# 		print group2	
		# 		group1 = re.match(r".*(\[.*\])(\*)$","".join(regex_list)).group(1)	
		# 	for j in group2:
		# 		regex_list.pop()
		# 	for k in group1:
		# 		regex_list.pop()
		#  	self.reverse_list += group1+group2 
		#  	if len(regex_list)> 1:
		#  		# ipdb.set_trace()
		#  		# regex_list = self.remove_pattern(r"\{.*\}", regex_list)
		# 		regex_seq_finder.regex_reverse(self, "".join(regex_list), nuctype="DNA")

		return self.reverse_list	

	def remove_pattern(pattern, regexl):
		if re.match(r".*(\(.*\))(pattern)$","".join(regexl)):
			# ipdb.set_trace()
			group2 = re.match(r".*(\(.*\))(pattern)$","".join(regexl)).group(2)
			print group2	
			group1 = re.match(r".*(\(.*\))(pattern)$","".join(regexl)).group(1)

		if re.match(r".*([A-Z]{1})(\{.*\})$","".join(regexl)):
			group2 = re.match(r".*([A-Z]{1})(pattern)$","".join(regexl)).group(2)
			print group2	
			group1 = re.match(r".*([A-Z]{1})(pattern)$","".join(regexl)).group(1)

		if re.match(r".*(\[.*\])(\{.*\})$","".join(regexl)):
			group2 = re.match(r".*(\[.*\])(pattern)$","".join(regexl)).group(2)
			print group2	
			group1 = re.match(r".*(\[.*\])(pattern)$","".join(regexl)).group(1)	
		for j in group2:
			regexl.pop()
		for k in group1:
			regexl.pop()
	 	self.reverse_list += group1+group2
	 	return regexl



	def find_subseq(self, sequence, regex, number_of_match, position_of_match, nuctype="DNA", overlap=False):
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
				# if nuctype == 'RNA':
				# 	value = IUPACData.ambiguous_rna_values[nt]
				# 	if len(value) == 1: 
				# 		pattern += value 
				# 	else: 
				# 		pattern += '[%s]' % value
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
	print("test")
	print(regex_seq_finder(object).regex_complement(r"AW{1,11}(CG){2,22}"))
	# print(regex_seq_finder(object).find_subseq("ATCTTTTTATTTCGCGCGGGGAAA",r"AW{1,10}(CG){1,10}", True, False ))
	val=regex_seq_finder(object).regex_reverse(r"AW{1,11}ACTTSS(CG){2,22}ACCTTA(CT))")
	print val 
	print(regex_seq_finder(object).regex_complement(val))