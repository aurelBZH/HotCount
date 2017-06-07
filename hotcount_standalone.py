#/usr/bin/python
#coding: utf-8
from __future__ import print_function
import argparse
import ConfigParser
import scipy
from hotcount import *
import scipy.stats as stats
import logging
import graypy
from logging.handlers import RotatingFileHandler

version = 1.0


class stand_alone(object):
	"""docstring for HotCount"""
	def __init__(self):
		# self.file_path = path
		# self.design_file = design_file
		# self.analysis_type = analysis_type
		# self.filetype = filetype
		self.analyse_result = ""
		self.stat_result = ""
		# self.pvalue = pvalue
		# self.mutation= mutation.split(",")
		# if output_file !="default" :
		# 	self.output_file =output_file

	# def choose_analysis(self):
	# 	if self.output_file:
	# 		try:
	# 			op_file = open(self.output_file, "w")
	# 		except Exception as e:
	# 			raise Exception("there is a problem with the file opening")
			
			
	# 		if self.analysis_type == "ALL":
	# 			self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
	# 			op_file.write(str(self.analyse_result))
	# 			self.stat_result = global_stat(analyse_result, self.pvalue,self.mutation)
	# 			op_file.write(str(self.stat_result))

	# 		elif self.analysis_type == "count":
	# 			self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
	# 			op_file.write(str(self.analyse_result))
	# 		elif self.analysis_type == "stat":
	# 			self.stat_result =global_stat(self.analyse_result, self.pvalue, self.mutation)
	# 			op_file.write(str(self.stat_result))
	# 		else:
	# 			raise ValueError("you have entered a false param. Param can only be either ALL, analysis or stat")	

	# 		op_file.close()	
	# 	else:

	# 		if self.analysis_type == "ALL":
	# 			self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
	# 			print(self.analyse_result)

	# 			self.stat_result = global_stat(self.analyse_result, self.pvalue, self.mutation)
	# 			print(self.stat_result)

	# 		elif self.analysis_type == "count":
	# 			self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
	# 			print(self.analyse_result)
	# 		elif self.analysis_type == "stat":
	# 			self.stat_result = global_stat(self.analyse_result, self.pvalue, self.mutation)
	# 			print(self.stat_result)
	# 		else:
				# raise ValueError("you have entered a false param. Param can only be either ALL, analysis or stat")	
		
	def count(self, design_file, path, filetype, output_file="default"):
		print (output_file)
		if output_file!=None:
			print("output_file")
			try:
				op_file = open(output_file, "w")
			except Exception as e:
				raise Exception("there is a problem with the file opening")	

			self.analyse_result = analysis(filetype, path, design_file)
			op_file.write(str(self.analyse_result))				
			op_file.close()	
		else:
			self.analyse_result = analysis(filetype, path, design_file)
			print(self.analyse_result)


	def ALL(self,design_file, path, filetype, pvalue, mutation, sample, output_file="default", controle):
		if output_file!=None:
			try:
				op_file = open(output_file, "w")
			except Exception as e:
				raise Exception("there is a problem with the file opening")
			self.analyse_result = analysis(filetype, file_path, design_file)
			op_file.write(str(self.analyse_result))
			self.stat_result = global_stat(self.analyse_result, pvalue,mutation, sample)
			op_file.write(str(to_csv(self.stat_result)))
			op_file.close()	
		else:
			self.analyse_result = analysis(filetype, file_path, design_file)
			print(self.analyse_result)

			self.stat_result = global_stat(self.analyse_result, pvalue, mutation)
			print(self.stat_result)


	def stat(self, countfile, pvalue, sample, mutation, output_file="default"):
		analyse_result = to_dict(countfile)
		if output_file!=None:
			try:
				op_file = open(self.output_file, "w")
			except Exception as e:
				raise Exception("there is a problem with the file opening")
			self.stat_result =global_stat(analyse_result, pvalue, mutation)
			op_file.write(str(to_csv(self.stat_result)))
			op_file.close()	
		else:
			self.stat_result = global_stat(analyse_result, pvalue, mutation)
			print(to_csv(self.stat_result))

def analysis(tmp_filetype, tmp_path, tmp_design_file):

	if tmp_filetype == "FASTQ":
		FQanalyse = analysisFQ()
		FQanalyse.get_file(tmp_path)
		analyse_result = FQanalyse.count_read(tmp_design_file)
		print(analyse_result)
	elif tmp_filetype == "BAM":
		BAManalyse = analysisBAM()
		BAManalyse.get_file(self.file_path)
		analyse_result = BAManalyse.count_read(tmp_design_file)
	return analyse_result	


def global_stat(analyse_result, pvalue,mutation):
	statistic = statistics(analyse_result, pvalue, *mutation)
	statistic.create_contingency_table()
	fisher = statistic.apply_fisher_test()
	results = statistic.check_positiv_sample()
	return results		

def to_csv(dicti, output_file ="default"):
	cpt=0
	if output_file == "default":
		output_file = "result.txt"
	with open(output_file, 'wb') as f:
		for i,j in dicti.iteritems(): 
			j["sample"]=i
			w = csv.DictWriter(f,j)
			if cpt == 0 :
			    w.writeheader()
			w.writerow(j)
			cpt+=1

def to_dict(file):
	with open(file, 'rb') as f:
		val = csv.DictReader(f)
		for row in val:
			print (row)
    	return val


if __name__ == '__main__':

	config = ConfigParser.SafeConfigParser()
	parser = argparse.ArgumentParser()

	sub = parser.add_subparsers(help="type of analysis",dest="analysis_type")
	all_parser = sub.add_parser("ALL")
	all_parser.add_argument('--designfile',required=True,help='path to the design file, the design file is the file containing the variant')
	all_parser.add_argument('--path',required=True, help='where the sample file are stored')
	all_parser.add_argument('-f','--filetype', default='FASTQ', help='file type to process')
	all_parser.add_argument("-o","--output", help="result file")
	all_parser.add_argument("-p","--pvalue", default=0.001, help="pvalue", type=int)
	all_parser.add_argument("-s", "--sample", default=6, help="max number of positive samples")
	all_parser.add_argument("-m", "--mutation", required= True,  help=" mutation to be analysed, separated by a ',' " )	
	all_parser.add_argument("-a","--controle",default="all", help="controle regex")

	count_parser = sub.add_parser("count")
	count_parser.add_argument('--designfile',required=True,help='path to the design file, the design file is the file containing the variant')
	count_parser.add_argument('--path',required=True, help='where the sample file are stored')
	count_parser.add_argument('-f','--filetype', default='FASTQ', help='file type to process')
	count_parser.add_argument("-o","--output", help="result file")
	
	stat_parser = sub.add_parser("stat")
	stat_parser.add_argument("-c","--countfile", required=True, help="file with count result")
	stat_parser.add_argument("-p","--pvalue", default=0.001, help="pvalue", type=int)
	stat_parser.add_argument("-s", "--sample", default=6, help="max number of positive samples")
	stat_parser.add_argument("-m", "--mutation", required= True,  help=" mutation to be analysed, separated by a ',' " )
	stat_parser.add_argument("-o","--output", help="result file")
	stat_parser.add_argument("-a","--controle", default="all", help="controle regex")


	

	args = vars(parser.parse_args())

	if "designfile" in args:
		pass
		config.read(args["designfile"])
		print (args["designfile"])

	

	if "mutation" in args:
		mutation =args["mutation"]
	else:
		mutation ="default"
	print(args)

	ipdb.set_trace()

	hotcount_stda = stand_alone()

	if args["analysis_type"]=="ALL":
		hotcount_stda.ALL(config,args["path"],args["filetype"], args["pvalue"], mutation, args["sample"], args["output"], args["controle"])
	elif args["analysis_type"]=="count":
		hotcount_stda.count(config, args["path"], args["filetype"], args["output"])
	else :	
		hotcount_stda.stat(args["countfile"], args["pvalue"], args["sample"], args["mutation"], args["output"], args["controle"])
