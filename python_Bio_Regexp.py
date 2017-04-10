#/usr/bin/python
import re
from Bio.Data import IUPACData
import ipdb

class regex_seq_finder(object):
	"""docstring for regex_seq_finder"""
	def __init__(self):
		self.sequence = ""
		self.regex_subseq = ""
		self.nuctype = "DNA"
		self.find_subseq_result = []
		self.complement = ""
		self.reverse_str=""


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
			return self.reverse_str
		
		# for i in reversed(regex_list):
			# ipdb.set_trace()
		i=regex_subseq[-1]
		print "beforet"+regex_subseq

		print i
		if re.match(r"[A-Z]",i) or i == ".":
			print "test "+i
	 		self.reverse_str += i
	 		regex_subseq =regex_subseq[:-1]
 			regex_seq_finder.regex_reverse(self, regex_subseq, nuctype="DNA")
	 	print regex_subseq		
	 	print self.reverse_str
		if i == '}' :
			# todoreverse group1
			if re.match(r".*(\(.*?\))(\{.*?\})$",regex_subseq):
				group2 = re.match(r".*(\(.*?\))(\{.*?\})$",regex_subseq).group(2)
				group1 = re.match(r".*(\(.*?\))(\{.*?\})$",regex_subseq).group(1)
				print group1
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
		 	self.reverse_str += group1+group2 
		 	regex_seq_finder.regex_reverse(self, regex_subseq, nuctype="DNA")

	 	if i == ')' :
	 		cpt=0
	 		tmp_regex=""
	 		list_regex = list(regex_subseq)
	 		for j in reversed(list_regex):
		 		ipdb.set_trace()
	 			if j ==')':
	 				cpt+=1
	 			if j == '(':
	 				cpt -=1
	 			if cpt == 0:
					regex_seq_finder.regex_reverse(self, tmp_regex, nuctype="DNA")
	 				break
	 			tmp_regex+j
	 		ipdb.set_trace()	
			self.reverse_str += "("
			tmp_regex_grp = re.escape(tmp_regex)+'$'
			regex_subseq = re.sub(tmp_regex_grp,'',regex_subseq)
			# regex_seq_finder.regex_reverse(self, regex_subseq, nuctype="DNA")
		if i == '(':
			self.reverse_str += ")"
		if i == ']' :	
			if re.match(r".*(\[.*?\])$", regex_subseq):
				group = re.match(r".*(\[.*?\])",regex_subseq).group(1)
				self.reverse_str += group
				tmp_regex_grp = re.escape(group)+'$'
				regex_subseq = re.sub(tmp_regex_grp,'',regex_subseq)
				regex_seq_finder.regex_reverse(self, regex_subseq, nuctype="DNA")

	
		if i == '*' :
			if re.match(r".*(\(.*?\))(\*|\+)$",regex_subseq):
				group2 = re.match(r".*(\(.*?\))(\*|\+)$",regex_subseq).group(2)
				group1 = re.match(r".*(\(.*?\))(\*)$",regex_subseq).group(1)

			elif re.match(r".*(\.|[A-Z]{1})(\*)$",regex_subseq):
				group2 = re.match(r".*(\.|[A-Z]{1})(\*|\+)$",regex_subseq).group(2)
				group1 = re.match(r".*(\.|[A-Z]{1})(\*|\+)$",regex_subseq).group(1)

			elif re.match(r".*(\[.*?\])(\*)$",regex_subseq):
				group2 = re.match(r".*(\[.*?\])(\*|\+|\+)$",regex_subseq).group(2)
				group1 = re.match(r".*(\[.*\])(\*|\+)$",regex_subseq).group(1)	

			tmp_regex_grp2 = re.escape(group2)+'$'
			tmp_regex_grp1 = re.escape(group1)+'$'
			regex_subseq = re.sub(tmp_regex_grp2,'',regex_subseq)
			regex_subseq = re.sub(tmp_regex_grp1,'',regex_subseq)
		 	self.reverse_str += group1+group2 

			regex_seq_finder.regex_reverse(self, regex_subseq, nuctype="DNA")


		return self.reverse_str	

	# def remove_pattern(pattern, regexl):
	# 	if re.match(r".*(\(.*\))(pattern)$","".join(regexl)):
	# 		# ipdb.set_trace()
	# 		group2 = re.match(r".*(\(.*\))(pattern)$","".join(regexl)).group(2)
	# 		# print group2	
	# 		group1 = re.match(r".*(\(.*\))(pattern)$","".join(regexl)).group(1)

	# 	if re.match(r".*([A-Z]{1})(\{.*\})$","".join(regexl)):
	# 		group2 = re.match(r".*([A-Z]{1})(pattern)$","".join(regexl)).group(2)
	# 		# print group2	
	# 		group1 = re.match(r".*([A-Z]{1})(pattern)$","".join(regexl)).group(1)

	# 	if re.match(r".*(\[.*\])(\{.*\})$","".join(regexl)):
	# 		group2 = re.match(r".*(\[.*\])(pattern)$","".join(regexl)).group(2)
	# 		# print group2	
	# 		group1 = re.match(r".*(\[.*\])(pattern)$","".join(regexl)).group(1)	
	# 	for j in group2:
	# 		regexl.pop()
	# 	for k in group1:
	# 		regexl.pop()
	#  	self.reverse_str += group1+group2
	#  	return self.reverse_str



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
				if nuctype == 'RNA':
					value = IUPACData.ambiguous_rna_values[nt]
					if len(value) == 1: 
						pattern += value 
					else: 
						pattern += '[%s]' % value
			elif re.match(r"[^A-Z]", nt):
				pattern += nt 
		self.find_subseq_result.append(pattern)
		if number_of_match == True :
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
# 	print("test")
# 	# print(regex_seq_finder(object).regex_complement(r"AW{1,11}(CG){2,22}"))
# 	# print(regex_seq_finder(object).find_subseq("ATCTTTTTATTTCGCGCGGGGAAA",r"AW{1,10}(CG){1,10}", True, False ))
	val=regex_seq_finder().regex_reverse(r"S(CG(TT))A")
	print val 
# 	# print(regex_seq_finder(object).regex_complement(val))