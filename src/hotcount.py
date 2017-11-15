#/usr/bin/env python 2.7
# coding: utf-8 
from __future__ import print_function

import gzip
from multiprocessing import Process, Queue, current_process

import glob2
import ipdb
import scipy.stats as stats
from Bio import SeqIO

from lib.python_bio_regexp.src import python_Bio_Regexp
from logsystem import *
from src.lib.pybam.src import pybam


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

    def count_read(self,design_dict,nuctype, IUPAC):
        """
        an abstract method
        implemented in the subclass*
        :param design_dict: the file with
        """
        raise NotImplementedError

    def create_regexp_dict(self, design_dict, nuctype, IUPAC):
        """
        a function to create a dictionnary of regexp with forwazr and reverse

        :param design_dict: dictionnary of regex to treat
        :param nuctype: type of nucleic acid
        :param IUPAC: use of iupac data
        :return: a dictionnary of regex with compressed regexp
        """
        logger.debug(design_dict)
        regexp_dict = {}
        for key , value in design_dict.iteritems():
            regexp_dict[key]={}
            regexp_dict[key]["forward"]=regex_seq_finder().create_pattern(value,nuctype,IUPAC)
            regexp_dict[key]["reverse"] = regex_seq_finder().create_pattern(regex_seq_finder().regex_reverse(value), nuctype, IUPAC)
        return regexp_dict

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

    def create_regexp_dict(self, design_dict, nuctype, IUPAC):
        """

        :param design_dict: dictionnary of regex
        :param nuctype: nucleic acid
        :param IUPAC: use of IUPAC code
        :return:
        """
        regexp_dict=super(AnalysisFQ, self).create_regexp_dict(design_dict, nuctype, IUPAC)
        return regexp_dict

    def count_read(self, design_dict,nuctype, IUPAC):
        """
        this function parrallelise the parsing with multiprocess
        :param design_dict: file containing the design regex
        :param nuctype : type of nucleic acid
        :param IUPAC : use IUPAC value or not
        """
        return_queue =Queue()
        compiled_regexp_dict = self.create_regexp_dict(design_dict, nuctype, IUPAC)

        procl = [Process(name=file, target=self.count_by_file, args=(compiled_regexp_dict, IUPAC, file, return_queue))for file in self.file_list]

        for proc in procl:
            proc.start()

        for proc in procl:
            proc.join()

            name, value = return_queue.get(True)
            self.analyse_results[name]= value

        return self.analyse_results


    def count_by_file(self, design_dict, IUPAC, file ,return_queue):
        """
        this method use pythonBioRegex library to count the number of match between
        design regex and sequence in file in FastQ
        :param design_dict: the hash of regexp to analyse
        :param IUPAC: use or don't use IUPAC code
        :param file: file name to treat
        :param return_queue: a queue to return results
        :return:
        """
        logger.info("treat %s" % file)
        self.analyse_results[file] = {}
        mutation_number_file_variant = {}
        if file.endswith(".gz"):
            fileContent = gzip.open(file, "rt")

        elif file.endswith(".gz") != True:
            fileContent = open(file, "rb")
        records = list(SeqIO.QualityIO.FastqGeneralIterator(fileContent))
        for name, design in design_dict.iteritems():
            mutation_number_by_var_val = 0
            mut_number = 0
            for record in records:
                tmp_mut = mut_number

                if (regex_seq_finder().find_subseq(str(record), design["forward"],True, IUPAC, False, False, True)[
                        1] and mut_number - 1 != tmp_mut) or (
                    regex_seq_finder().find_subseq(str(record), design["reverse"], True, IUPAC, False, False, True)[
                        1] and mut_number - 1 != tmp_mut):
                    mut_number = mut_number + 1
            mutation_number_by_var_val = mutation_number_by_var_val + mut_number
            mutation_number_file_variant[name] = mutation_number_by_var_val
        logger.debug("end of count %s"%file)

        return_queue.put((current_process().name,mutation_number_file_variant))

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

    def create_regexp_dict(self, design_dict, nuctype, IUPAC):
        regexp_dict = super(AnalysisBAM, self).create_regexp_dict(design_dict, nuctype, IUPAC)
        return regexp_dict


    def count_read(self, design_dict, IUPAC):
        """
        this method use pythonBioRegex library to count the number of match between
        design regex and sequence in file in BAM
        :param design_dict: file containing the design regex
        """
        read_set = set()
        for file in self.file_list:
            mutation_number_file_variant = self.count_by_file(design_dict, IUPAC,file, read_set)
            #samfile = pysam.AlignmentFile(file)
            #logger.info("treat %s"%file)
            #self.analyse_results[file] = {}
            #mutation_number_file_variant = {}
            #for name, design in design_dict.iteritems():
            #	mutation_number_by_var_val = 0
            #	reverse_design = regex_seq_finder().regex_reverse_complement(design)
            #	mut_number = 0
            #	for read in samfile.fetch(until_eof=True):
    #
                    #str_read=str(read).split("\t")

                    #if str_read[0] not in read_set:
                    #	logger.debug(re.match(design, str_read[9]))
                    #	logger.debug(design)
                    #	logger.debug(str_read[9])
                    #	ipdb.set_trace()
                    #	if regex_seq_finder().find_subseq(str(str_read[9]),design,IUPAC,False, False, True)[1]:
                    #		logger.debug(str_read[9])
                    #		ipdb.set_trace()

                    #		mut_number = mut_number+1
                    #	elif regex_seq_finder().find_subseq(str(str_read[9]),reverse_design, IUPAC,False, False, True)[1]:
                    #		mut_number = mut_number+1
                    #read_set.add(str_read[0])
                #mutation_number_by_var_val = mutation_number_by_var_val + mut_number
                #mutation_number_file_variant[name] = mutation_number_by_var_val
            self.analyse_results[file] = mutation_number_file_variant
        return self. analyse_results

    def count_by_file(self, design_dict, IUPAC, file, read_set):
        samfile = pysam.AlignmentFile(file)
        logger.info("treat %s" % file)
        self.analyse_results[file] = {}
        mutation_number_file_variant = {}
        for name, design in design_dict.iteritems():
            mutation_number_by_var_val = 0
            reverse_design = regex_seq_finder().regex_reverse_complement(design)
            mut_number = 0
            for read in samfile.fetch(until_eof=True):

                str_read = str(read).split("\t")

                if str_read[0] not in read_set:
                    logger.debug(re.match(design, str_read[9]))
                    logger.debug(design)
                    logger.debug(str_read[9])
                    if regex_seq_finder().find_subseq(str(str_read[9]), design, IUPAC, False, False, True)[1]:
                        logger.debug(str_read[9])

                        mut_number = mut_number + 1
                    elif regex_seq_finder().find_subseq(str(str_read[9]), reverse_design, IUPAC, False, False, True)[1]:
                        mut_number = mut_number + 1
                read_set.add(str_read[0])
            mutation_number_by_var_val = mutation_number_by_var_val + mut_number
            mutation_number_file_variant[name] = mutation_number_by_var_val

            return mutation_number_by_var_val


class statistics(object):
	"""docstring for statistics
	a class for statistics analysis

	"""
	def __init__(self, count_table, pvalue, sample, controle, depth, mutation):
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
		self.depth = depth
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
				if (self.controle == "all") and (int(self.count_table[sample][self.controle])>= int(sum(self.count_table[sample].values()))):
						raise ValueError("all is less than the sum of all the variant")
				if int(self.count_table[sample][mut_lower]) > int(self.count_table[sample][self.controle]):
					raise ValueError("Nombre de read porteur de l'expression régulière \"controle\" "+self.controle+ " infèrieur" +" ("+self.count_table[sample][self.controle]+")"+
						"au nombre de read muté"+" "+mut_lower+" ("+self.count_table[sample][mut_lower]+") "+". Rechercher une erreur dans l'expression régulière.")
				if int(self.count_table[sample]["all"])<self.depth:
					raise ValueError("the sequencing depth at the mutation position is insufficient. ")
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
