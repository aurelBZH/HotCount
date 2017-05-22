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
	def __init__(self, design_file, path, analysis_type, filetype, pvalue, mutation,  output_file="default"):
		self.file_path = path
		self.design_file = design_file
		self.analysis_type = analysis_type
		self.filetype = filetype
		self.analyse_result = ""
		self.stat_result = ""
		self.pvalue = pvalue
		self.mutation= mutation.split(",")
		if output_file !="default" :
			self.output_file =output_file

	def choose_analysis(self):
		if self.output_file:
			try:
				op_file = open(self.output_file, "w")
			except Exception as e:
				raise Exception("there is a problem with the file opening")
			
			
			if self.analysis_type == "ALL":
				self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
				op_file.write(str(self.analyse_result))
				self.stat_result = global_stat(analyse_result, self.pvalue,self.mutation)
				op_file.write(str(self.stat_result))

			elif self.analysis_type == "analysis":
				self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
				op_file.write(str(self.analyse_result))
			elif self.analysis_type == "stat":
				self.stat_result =global_stat(self.analyse_result, self.pvalue, self.mutation)
				op_file.write(str(self.stat_result))
			else:
				raise ValueError("you have entered a false param. Param can only be either ALL, analysis or stat")	

			op_file.close()	
		else:

			if self.analysis_type == "ALL":
				self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
				print(self.analyse_result)

				self.stat_result = global_stat(self.analyse_result, self.pvalue, self.mutation)
				print(self.stat_result)

			elif self.analysis_type == "analysis":
				self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
				print(self.analyse_result)
			elif self.analysis_type == "stat":
				self.stat_result = global_stat(self.analyse_result, self.pvalue, self.mutation)
				print(self.stat_result)
			else:
				raise ValueError("you have entered a false param. Param can only be either ALL, analysis or stat")	
		



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

def to_csv(dict):
	cpt=0
	with open('mycsvfile.csv', 'wb') as f:
        for i,j in count.iteritems(): 
            j["sample"]=i
            w = csv.DictWriter(f,j)
            if cpt == 0 :
    	        w.writeheader()
            w.writerow(j)
            cpt+=1

def to_dict():
	with open('mycsvfile.csv', 'rb') as f:
    ...:   val = csv.DictReader(f)
    ...:   for row in val:
    ...:         print row
    return val


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	config = ConfigParser.SafeConfigParser()

	parser.add_argument('--designfile',required=True,help='path to the design file, the design file is the file containing the variant')
	parser.add_argument('path', help='where the sample file are stored')
	parser.add_argument('-t','--analysistype', default='ALL', help='type of analysis to be processed')
	parser.add_argument('-f','--filetype', default='FASTQ', help='file type to process')
	parser.add_argument("-o","--output", help="result file")
	parser.add_argument("-p","--pvalue", default=0.001, help="pvalue", type=int)
	parser.add_argument("-s", "--sample", default=6, help="max number of positive samples")
	parser.add_argument("-m", "--mutation", required= True,  help=" mutation to be analysed, separated by a ',' " )	

	args = parser.parse_args()
	config.read(args.designfile)
	if args.mutation:
		mutation =args.mutation
	else:
		mutation ="default"
	hotcount_stda = stand_alone(config, args.path, args.analysistype, args.filetype, args.pvalue, mutation, args.output)
	hotcount_stda.choose_analysis()
