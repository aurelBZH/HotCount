#/usr/bin/python
# coding: utf-8 
from __future__ import print_function
from Bio import SeqIO
import glob2
import pysam
import scipy.stats as stats
from python_Bio_Regexp import *
from logsystem import *
import collections


class Analysis(object):
	"""
	super class for analysisFQ and analysisBAM
	"""

	def __init__(self):
		self.file_list=[]
		self.design_file = ""
		self.analyse_results = {}
		self.file_extentions = []
		self.file_list=[]

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
			raise Exception("the path is incorrect or the file extension"+ 
				"is bad,in fact the program can't find the file ")
		self.file_list = file_list                                                                            
		return self.file_list
		# pysam.AlignmentFile("ex1.bam", "rb")
        
	def count_read( design_file ): 
		"""
		an abstract method
		implemented in the subclass*
		:param design_file: the file with
		"""
		raise NotImplementedError


class AnalysisFQ(Analysis):
	"""
	docstring for analysisFQ
	a class inhériting from the analysis class allow to analyse FastQ file
	"""
	def __init__(self):
		super(AnalysisFQ, self).__init__()
		self.file_extentions = ['*.fq.*', '*.fastq','*.fq','*.FASTQ','*.fastq.gz', '.fq.gz', '.FASTQ.gz']
		self.file_list=[]
		self.analyse_results ={}


	def get_file(self, path):
		"""
		this method use the super class method
		:param path: path where the file are stored
		"""
		self.file_list = super(AnalysisFQ, self).get_file(path)


	def count_read(self, design_file):
		"""
		this method use pythonBioRegex library to count the number of match between
		design regex and sequence in file in FastQ
		:param design_file: file containing the design regex
		"""
		design_dict=dict(design_file.items('mutation_design'))

		for file in self.file_list:
			logger.info("treat %s"%file)

			if file.endswith(".gz") :
				pass
			self.analyse_results[file] = {}
			mutation_number_file_variant = {}
			with open(file, "rU") as handle:
				records = list(SeqIO.parse(handle, "fastq"))
				logger.info("treat %s"%file)
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

class AnalysisBAM(Analysis):
    """
    docstring for analysisBAM
    a class inhériting from the analysis class allow to analyse BAM

    """
    def __init__(self):
        super(AnalysisBAM, self).__init__()
        self.file_extentions = ['*.bam.*', '*.BAM.*', '*.bam','*.BAM']
        self.file_list=[]
        self.analyse_results = {}


		
	def get_file(self, path):
		"""
		this method use the super class method
		:param path: path where the file are stored
		"""

		self.file_list = super(AnalysisBAM, self).get_file(path)

    def count_read(self, design_file):
        """
        this method use pythonBioRegex library to count the number of match between
        design regex and sequence in file in FastQ
        :param design_file: file containing the design regex
        """
        design_dict=dict(design_file.items('mutation_design'))
        for file in self.file_list:
            samfile = pysam.AlignmentFile(file)
            logger.info("treat %s"%file)
            self.analyse_results[file] = {}
            mutation_number_file_variant = {}
            for name, design in design_dict.iteritems():
                mutation_number_by_var_val = 0
                reverse_design = regex_seq_finder().regex_reverse_complement(design)
                mut_number = 0
                for read in samfile.fetch():
                    if regex_seq_finder().find_subseq(read,design,False, False, True)[1]:
                        mut_number = mut_number+1
                    elif regex_seq_finder().find_subseq(read,reverse_design,False, False, True)[1]:
                        mut_number = mut_number+1
                mutation_number_by_var_val= mutation_number_by_var_val + mut_number
                mutation_number_file_variant[name] = mutation_number_by_var_val
            self.analyse_results[file] = mutation_number_by_var_val
        return self. analyse_results

class statistics(object):
	"""docstring for statistics
	a class for statistics analysis

	"""
	def __init__(self, count_table, pvalue, sample, controle, *mutation):
		"""
		:param count_table: dictionary of count result
		:param pvalue: the pvalue to use for stat analysis
		:param samle: the number of positiv sample
		:param controle: regex to use as controle during the stat analysis
		:param mutation: mutation to analyse
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
		"""
		a methode to create a contingency table
		usable for fisher test.

		"""
		for sample in self.count_table:
			# ipdb.set_trace()
			logger.info(self.mutation)
			for j in self.mutation:
				try:
					if self.count_table[sample][j] > self.count_table[sample][self.controle] :
						raise ValueError("Nombre de read porteur de l'expression régulière \"controle\" infèrieur" +
							"au nombre de read muté. Rechercher une erreur dans l'expression régulière.")
					self.contingency_table[sample] = [self.count_table[sample][j],self.count_table[sample][self.controle]]

				except ValueError:
					raise ValueError("problem in value error ")
					# print(self.contingency_table)		
		return self.contingency_table 
			
	def apply_fisher_test(self):
		"""
		a method to apply fisher test on. It's a statistic test to check if a test
		"""
		fisher_hash_result={}
		# ipdb.set_trace()
		for sample1, count_value1  in self.contingency_table.iteritems():
			fisher_hash_result[sample1] = {}
			for sample2, count_value2 in self.contingency_table.iteritems():
				if sample1 != sample2:

					fisher_result = stats.fisher_exact([count_value1, count_value2],"greater")
					fisher_hash_result[sample1][sample2] = fisher_result[1]
		self.fisher_matrix = fisher_hash_result
		return self.fisher_matrix	

	def check_positiv_sample(self):
		"""
		a method to check sample positivity
		"""
		positivity = False
		pvalue_tab= []
		sample_pvalue_dict = {}
		for sample_name, dict_result in self.fisher_matrix.iteritems() :
			sample_pvalue_dict = {sample_name: []}
			for sample2, pvalue in dict_result.iteritems():
				sample_pvalue_dict[sample_name].append(pvalue)
			pvalue_tab.append(sample_pvalue_dict)
		significativ_value_dict = collections.OrderedDict()
		for sample_res in pvalue_tab:
			nb_neg=0
			for sample_name, pvalues in sample_res.iteritems():
				minimal_pvalue=self.pvalue
				minimal_samplename =""
				for value in pvalues:
					if value <minimal_pvalue:
						minimal_pvalue = value
						minimal_samplename = sample_name
				if minimal_pvalue<self.pvalue:
					significativ_value_dict[minimal_samplename] = minimal_pvalue
		return significativ_value_dict


if __name__ == '__main__':
	AnalysisFQ().get_file("resources/")
	AnalysisFQ.count_read()
	statistics()