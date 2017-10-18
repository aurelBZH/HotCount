#/usr/bin/env python 2.7
# coding: utf-8 
from __future__ import print_function
from Bio import SeqIO
import glob2
import pysam
import scipy.stats as stats
from python_Bio_Regexp import *
from logsystem import *

class Analysis(object):
	"""
	super class for analysisFQ and analysisBAM
	"""

	def __init__(self):
		self.file_list=[]
		self.design_dict = ""
		self.analyse_results = {}
		self.file_extentions = []

	def get_file(self, path):
		"""
		a method to get the file list from the path 
		:param path: path where the file are stored
		"""
		for extention in self.file_extentions :
			file_list = glob2.glob(path+extention)
			if len(file_list) != 0 :
					break
		if len(file_list) == 0 :
			raise Exception("the path is incorrect or the file extension "+
		 		"is bad,in fact the program can't find the file ")

		self.file_list = file_list
		return file_list
		# pysam.AlignmentFile("ex1.bam", "rb")
        
	def count_read(design_dict):
		"""
		an abstract method
		implemented in the subclass*
		:param design_dict: the file with
		"""
		raise NotImplementedError


class AnalysisFQ(Analysis):
	"""
	docstring for analysisFQ
	a class inhériting from the analysis class allow to analyse FastQ file
	"""
	def __init__(self):
		super(AnalysisFQ, self).__init__()
		self.file_extentions = ['*.fq.*', '*.fastq','*.fq','*.FASTQ','*.fastq.gz', '.fq.gz', '.FASTQ.gz', ".FQ", ".FQ.gz"]
		self.file_list=[]
		self.analyse_results ={}


	def get_file(self, path):
		"""
		this method use the super class method
		:param path: path where the file are stored
		"""
		self.file_list = super(AnalysisFQ, self).get_file(path)


	def count_read(self, design_dict):
		"""
		this method use pythonBioRegex library to count the number of match between
		design regex and sequence in file in FastQ
		:param design_dict: file containing the design regex
		"""
		for file in self.file_list:
			logger.info("treat %s"%file)

			if file.endswith(".gz") :
				pass
			self.analyse_results[file] = {}
			mutation_number_file_variant = {}
			with open(file, "rU") as handle:
				records = list(SeqIO.parse(handle, "fastq"))
				for name, design in design_dict.iteritems():
					mutation_number_by_var_val = 0
					reverse_design = regex_seq_finder().regex_reverse_complement(design)
					mut_number = 0
					for record in records :
						tmp_mut = mut_number
						if (regex_seq_finder().find_subseq(str(record.seq),design,False, False, True)[1] and mut_number-1!=tmp_mut) or (regex_seq_finder().find_subseq(str(record.seq),reverse_design,False, False, True)[1] and mut_number-1!=tmp_mut):
							mut_number = mut_number+1
#						elif regex_seq_finder().find_subseq(str(record.seq),reverse_design,False, False, True)[1] & mut_number-1==tmp_mut:
#							mut_number = mut_number+1
					mutation_number_by_var_val= mutation_number_by_var_val + mut_number
					mutation_number_file_variant[name] = mutation_number_by_var_val
			self.analyse_results[file] = mutation_number_file_variant
		return self. analyse_results

class AnalysisBAM(Analysis):
	"""
	docstring for analysisBAM
	a class inhériting from the analysis class allow to analyse BAM

	"""
	def __init__(self):
		super(AnalysisBAM, self).__init__()
		self.file_extentions = ['*.bam','*.BAM']
		self.file_list=[]
		self.analyse_results = {}


		
	def get_file(self, path ):
		"""
		this method use the super class method
		:param path: path where the file are stored
		"""

		self.file_list = super(AnalysisBAM, self).get_file(path)

	def count_read(self, design_dict):
		"""
		this method use pythonBioRegex library to count the number of match between
		design regex and sequence in file in BAM
		:param design_dict: file containing the design regex
		"""
		read_set = set()
		for file in self.file_list:
			samfile = pysam.AlignmentFile(file)
			logger.info("treat %s"%file)
			self.analyse_results[file] = {}
			mutation_number_file_variant = {}
			for name, design in design_dict.iteritems():
				mutation_number_by_var_val = 0
				reverse_design = regex_seq_finder().regex_reverse_complement(design)
				mut_number = 0
				for read in samfile.fetch(until_eof=True):

					str_read=str(read).split("\t")
					if str_read[0] not in read_set:
						if regex_seq_finder().find_subseq(str_read[9],design,False, False, True)[1]:
							mut_number = mut_number+1
						elif regex_seq_finder().find_subseq(str_read[9],reverse_design,False, False, True)[1]:
							mut_number = mut_number+1
					read_set.add(str_read[0])
				mutation_number_by_var_val = mutation_number_by_var_val + mut_number
				mutation_number_file_variant[name] = mutation_number_by_var_val
			self.analyse_results[file] = mutation_number_file_variant
		return self. analyse_results

class statistics(object):
	"""docstring for statistics
	a class for statistics analysis

	"""
	def __init__(self, count_table, pvalue, sample, controle, mutation):
		"""
		:param count_table: dictionary of count result
		:param pvalue: the pvalue to use for stat analysis
		:param sample: the number of positiv sample
		:param controle: regex to use as controle during the stat analysis
		:param mutation: mutation to analysecount_value2
		:type count_table: dictionary
		:type pvalue: float
		:type sample: int
		:type controle: string

		:type mutation: string
		"""
		self.contingency_table = {}
		self.count_table = count_table
		self.fisher_matrix =[]
		self.mutation = mutation
		self.pvalue =pvalue
		self.sample = sample
		self.controle = controle
		logger.info("begin stat")


	def create_contingency_table(self):
		"""count_value2
		a methode to create a contingency table
		usable for fisher test.

		"""
		for mut in self.mutation:
			mut_lower = mut.lower()
			if mut_lower not in self.contingency_table:
				self.contingency_table[mut_lower] = {}
			for sample in self.count_table:
				logger.info(" %s"%mut)

				if int(self.count_table[sample][mut_lower]) > int(self.count_table[sample][self.controle]) :
					if self.controle == "all" & int(self.count_table[sample][self.controle])>= sum(self.count_table[sample]):
						raise ValueError("all is less than the sum of all the variant")
					raise ValueError("Nombre de read porteur de l'expression régulière \"controle\" "+self.controle+ " infèrieur" +" ("+self.count_table[sample][self.controle]+")"+
						"au nombre de read muté"+" "+mut_lower+" ("+self.count_table[sample][mut_lower]+") "+". Rechercher une erreur dans l'expression régulière.")

				self.contingency_table[mut_lower][sample] = [self.count_table[sample][mut_lower],self.count_table[sample][self.controle]]
	#si controle =all (soustraire
		return self.contingency_table
			
	def apply_fisher_test(self):
		"""
		a method to apply fisher test on. It's a statistic test to check if a value is superior to the noise.
		each sample count is compared with all others sample count.

		"""
		fisher_hash_result={}
		# ipdb.set_trace()

		for mutation, contingency_table in self.contingency_table.iteritems():
			if mutation not in fisher_hash_result:
				fisher_hash_result[mutation]={}
				logger.debug(mutation)
			for sample1, count_value1  in contingency_table.iteritems():
				fisher_hash_result[mutation][sample1] = {}
				for sample2, count_value2 in contingency_table.iteritems():
					if sample1 != sample2:
						fisher_result = stats.fisher_exact([count_value1, count_value2],"greater")
						fisher_hash_result[mutation][sample1][sample2] = fisher_result[1]
		self.fisher_matrix = fisher_hash_result
		logger.info(self.fisher_matrix)
		return self.fisher_matrix

	def check_positiv_sample(self):
		"""
		a method to check sample positivity, check for each couple of
		sample if the pvalue is inferior to the chosen pvalue

		"""
		positivity = False
		sample_pvalue_dict = {}
		for mutation, dict_result in self.fisher_matrix.iteritems() :
			if mutation not in sample_pvalue_dict:
				sample_pvalue_dict[mutation]= {}
			for sample_name, result in dict_result.iteritems():
				sample_pvalue_dict[mutation][sample_name]=[]
				for sample2, pvalue in result.iteritems():
					sample_pvalue_dict[mutation][sample_name].append(pvalue)
		significativ_value_dict = {}
		for mutation, final_pvalue_dict in sample_pvalue_dict.iteritems():
			result_by_mutation=[]

			for sample_name, pvalues in final_pvalue_dict.iteritems():
				sorted_pvalue_list = sorted(pvalues, reverse=True)
				result_by_mutation.append([sample_name,sorted_pvalue_list])

			significativ_value_dict[mutation]=result_by_mutation


		return significativ_value_dict
