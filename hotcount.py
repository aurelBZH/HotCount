#/usr/bin/python
from Bio import SeqIO
import glob2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqUtils
import ipdb
import pysam
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
			if len(file_list) != 0 :
					break
		if len(file_list) == 0 :	 
			raise Exception("the path is incorrect or the file extension"+ 
				"is bad,in fact the program can't find the file ")
		self.file_list = file_list                                                                            

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

	# def statistics(self):
	# 	print("why")
	# 	print("test")
	
class analysisFQ(analysis):
	"""docstring for analysisFQ"""
	def __init__(self):
		super(analysisFQ, self).__init__()
		self.file_extentions = ['*.fastq.*', '*.fq.*', '*.fastq','*.fq','*.FASTQ']
		self.file_list=[]
		self.design_file = ""


		
	def get_file(self, path):
		super(analysisFQ, self).get_file(path)
		

	def count_read(self, design_file):
				design_dict={}
		for design_values in design_file.items('mutation_design'):
			design_dict[design_values[0]] = Seq(design_values[1]) 
           
		for file in self.file_list:
			self.analyse_results[file] = {}
			mutation_number_file_variant = []
			with open(file, "rU") as handle:
				records = list(SeqIO.parse(handle, filetype.lower()))
  
				for name, design in design_dict.iteritems():
					mutation_number_by_var_val = 0       
					print design
					reverse_design = python_Bio_Regexp().regex_reverse_complement(design)
					mut_number = 0
					for record in records : 
						if python_Bio_Regexp().find_subseq(record,design,False, False, True)[1]:
							mut_number = mut_number+1
						elif python_Bio_Regexp().find_subseq(record,reverse_design,False, False, True)[1]: 
							mut_number = mut_number+1 
					mutation_number_by_var_val= mutation_number_by_var_val + mut_number
					#ipdb.set_trace()
					mutation_number_file_variant.append([name,mutation_number_by_var_val])
			self.analyse_results[file] = mutation_number_file_variant 
			return self. analyse_results    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
	
		

class analysisBAM(analysis):
	"""docstring for analysisBAM"""
	def __init__(self):
		super(analysisBAM, self).__init__()
		self.file_extentions = ['*.bam.*', '*.BAM.*', '*.bam','*.BAM']
		self.file_list=[]
		self.design_file = ""


		
	def get_file(self, path):
		super(analysisBAM,self).get_file(path)

	def count_read(self, design_file):
		design_dict={}
		for design_values in design_file.items('mutation_design'):
			design_dict[design_values[0]] = Seq(design_values[1]) 

		for file in self.file_list:
			samfile = pysam.AlignmentFile(file)
			self.analyse_results[file] = {}
			mutation_number_file_variant = []

			for name, design in design_dict.iteritems():
				mutation_number_by_var_val = 0       
				print design
				mut_number = 0
				for read in samfile.fetch():
					if python_Bio_Regexp().find_subseq(record,design,False, False, True)[1]:
							mut_number = mut_number+1
					elif python_Bio_Regexp().find_subseq(record,reverse_design,False, False, True)[1]: 
							mut_number = mut_number+1 
				mutation_number_by_var_val= mutation_number_by_var_val + mut_number
					#ipdb.set_trace()
				mutation_number_file_variant.append([name,mutation_number_by_var_val])
			self.analyse_results[file] = mutation_number_file_variant 
			return self. analyse_results    


if __name__ == '__main__':
	analysisFQ().get_file("resources/")
	analysisFQ.count_read()
