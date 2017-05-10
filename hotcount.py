#/usr/bin/python
# coding: utf-8 
from __future__ import print_function
from Bio import SeqIO
import glob2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqUtils
import ipdb
import pysam
import scipy.stats as stats
from python_Bio_Regexp import *


class analysis(object):

	def __init__(self):
		self.file_list=[]
		self.design_file = ""
		self.analyse_results = {}
		self.file_extentions = []
		self.file_list=[]

	def get_file(self, path):
		print self.file_extentions
		for extention in self.file_extentions : 
			file_list = glob2.glob(path+extention)
			# print file_list
			if len(file_list) != 0 :
					break
		if len(file_list) == 0 :	 
			raise Exception("the path is incorrect or the file extension"+ 
				"is bad,in fact the program can't find the file ")
		self.file_list = file_list                                                                            
		return self.file_list
		# pysam.AlignmentFile("ex1.bam", "rb")
        
	def count_read( design_file ): 
		raise NotImplementedError
		# design_dict={}
		# for design_values in design_file.items('mutation_design'):
		# 	design_dict[design_values[0]] = Seq(design_values[1]) 
           
		# for file in self.file_list:
		# 	self.analyse_results[file] = {}
		# 	mutation_number_file_variant = []
		# 	with open(file, "rU") as handle:
		# 		records = list(SeqIO.parse(handle, filetype.lower()))
  
		# 		for name, design in design_dict.iteritems():
		# 			mutation_number_by_var_val = 0       
		# 			print design
		# 			mut_number = 0
		# 			for record in records : 
		# 				if len(SeqUtils.nt_search(str(record.seq),str(design)))>1:
		# 					mut_number = mut_number+1
		# 				elif len(SeqUtils.nt_search(str(record.seq),str(design.reverse_complement())))>1: 
		# 					mut_number = mut_number+1 
		# 			mutation_number_by_var_val= mutation_number_by_var_val + mut_number
		# 			#ipdb.set_trace()
		# 			mutation_number_file_variant.append([name,mutation_number_by_var_val])
		# 	self.analyse_results[file] = mutation_number_file_variant                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
	
					
			# elif filetype == "BAM":
			# 	pass
			# else : 
				# raise TypeError(" the file type you have chosen isn't processed by the program")
	
class analysisFQ(analysis):
	"""docstring for analysisFQ"""
	def __init__(self):
		super(analysisFQ, self).__init__()
		self.file_extentions = ['*.fastq.*', '*.fq.*', '*.fastq','*.fq','*.FASTQ']
		self.file_list=[]
		self.analyse_results ={}

		
	def get_file(self, path):
		self.file_list = super(analysisFQ, self).get_file(path)
		print(self.file_list)
		

	def count_read(self, design_file):
		design_dict={}
		for design_values in design_file.items('mutation_design'):
			design_dict[design_values[0]] = design_values[1] 	 

		for file in self.file_list:
			if file.endswith(".gz") :
				pass
			self.analyse_results[file] = {}
			mutation_number_file_variant = {}
			with open(file, "rU") as handle:
				records = list(SeqIO.parse(handle, "fastq"))
  
				for name, design in design_dict.iteritems():
					mutation_number_by_var_val = 0       
					print(design)
					reverse_design = regex_seq_finder().regex_reverse_complement(design)
					mut_number = 0
					for record in records : 
						if regex_seq_finder().find_subseq(str(record.seq),design,False, False, True)[1]:
							mut_number = mut_number+1
						elif regex_seq_finder().find_subseq(str(record.seq),reverse_design,False, False, True)[1]: 
							mut_number = mut_number+1 
					mutation_number_by_var_val= mutation_number_by_var_val + mut_number
					mutation_number_file_variant[name] = mutation_number_by_var_val
			self.analyse_results[file] = mutation_number_file_variant 
		return self. analyse_results

class analysisBAM(analysis):
	"""docstring for analysisBAM"""
	def __init__(self):
		super(analysisBAM, self).__init__()
		self.file_extentions = ['*.bam.*', '*.BAM.*', '*.bam','*.BAM']
		self.file_list=[]
		self.analyse_results = {}


		
	def get_file(self, path):
		self.file_list = super(analysisFQ, self).get_file(path)

	def count_read(self, design_file):
		design_dict=dict(design_file.items('mutation_design'))
		for design_values in design_dict:
			design_dict[design_values[0]] = Seq(design_values[1]) 
		for file in self.file_list:
			samfile = pysam.AlignmentFile(file)
			self.analyse_results[file] = {}
			mutation_number_file_variant = {}

			for name, design in design_dict.iteritems():
				mutation_number_by_var_val = 0       
				mut_number = 0
				for read in samfile.fetch():
					if regex_seq_finder().find_subseq(read,design,False, False, True)[1]:
							mut_number = mut_number+1
					elif regex_seq_finder().find_subseq(read,reverse_design,False, False, True)[1]: 
							mut_number = mut_number+1 
				mutation_number_by_var_val= mutation_number_by_var_val + mut_number
					#ipdb.set_trace()
				mutation_number_file_variant[name] = mutation_number_by_var_val
			self.analyse_results[file] = mutation_number_by_var_val 
			return self. analyse_results    

class statistics(object):
	"""docstring for statistics"""
	def __init__(self, count_table, *mutation):
		super(statistics, self).__init__()
		self.contingency_table = {}
		self.count_table = count_table
		self.fisher_matrix =[]
		self.mutation = mutation

	def create_contingency_table(self):
		for sample in self.count_table:
			# ipdb.set_trace()
			for j in self.mutation:
				try:
					self.contingency_table[sample] = [self.count_table[sample][j],self.count_table[sample]["all"]]

				except Exception as Valueerror:
					raise e
		print(self.contingency_table)		
		return self.contingency_table 
			
	def apply_fisher_test(self):
		fisher_hash_result={}
		# ipdb.set_trace()
		for sample1, count_value1  in self.contingency_table.iteritems():
			
			fisher_hash_result["sample1"] = {}
			for sample2, count_value2 in self.contingency_table.iteritems():
				if sample1 != sample2:
					fisher_result = stats.fisher_exact([count_value1, count_value2],"greater")
					fisher_hash_result[sample1][sample2] = fisher_result[1]
		print(fisher_result)
		self.fisher_matrix = fisher_hash_result
		return self.fisher_matrix	

	def check_positiv_sample(self):
		positivity = False
		pvalue_tab= []
		print (self.fisher_matrix)
		for sample_name, dict_result in self.fisher_matrix.iteritems() :
			sample_pvalue_dict = {sample_name: []}
			for sample2, pvalue in dict_result.iteritems():
				sample_pvalue_dict[sample_name].append(pvalue)
				pvalue_tab.append(sample_pvalue_dict)
		nb_neg=0		
		
		for sample_res in pvalue_tab:
			print (pvalue_tab)
			for sample_name, pvalue in sample_res.iteritems():
				if pvalue > 0.01:
					nb_neg +=1
			
			return "{0}is positiv with only {1} ".format(sample_name,nb_neg)					

if __name__ == '__main__':
	analysisFQ().get_file("resources/")
	analysisFQ.count_read()
	statistics()