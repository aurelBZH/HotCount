#/usr/bin/python
import re

class regex_seq_finder(object):
	"""docstring for regex_seq_finder"""
	def __init__(self, arg):
		self.arg = arg
		self.sequence = ""
		self.regex_subseq = ""
		self.type = "DNA"
		self.find_subseq_result = []
		

	
	def find_reverse_complement(regex_subseq, type):
		self.regex_subseq = regex_subseq


	def find_subseq(sequence, overlap=True, number_of_match, position_of_match):
			self.sequence = sequence	
			def nt_search(seq, subseq): 

	    pattern = '' 
    	for nt in subseq:
    		if re.match(r"[A-Z]",nt): 
	        	value = IUPACData.ambiguous_dna_values[nt] 
	        	if len(value) == 1: 
            		pattern += value 
          		else: 
	            	pattern += '[%s]' % value 
   
        self.find_subseq_result.append(pattern)
        if number_of_match ==True :
        	number_subseq = 0 
        	number_subseq+=1 if len(re.findall(pattern, self.sequence,overlap))>0 
    		self.find_subseq_result.append(number_subseq) 
	    	return self.find_subseq_result 
	    if position_of_match == True:
	    	matches=re.finditer(pattern, self.sequence)
	    	for match in matches :
	    		self.find_subseq_result.append(match.start)
	    	return self.find