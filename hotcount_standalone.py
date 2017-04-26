#/usr/bin/python

import argparse
import ConfigParser
import numpy
from hotcount import *

version = 1.0


class stand_alone(object):
	"""docstring for HotCount"""
	def __init__(self, design_file, path, analysis_type, filetype):
		self.file_path = path
		self.design_file = design_file
		self.analysis_type = analysis_type
		self.filetype = filetype
		self.analyse_result = ""
		self.stat_result = ""

	def choose_analysis(self):
		if self.analysis_type == "ALL":
			self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
			self.stat_result = statistics()

		elif self.analysis_type == "analysis":
			self.analyse_result = analysis(self.filetype, self.file_path, self.design_file)
		elif self.analysis_type == "stat":
			self.stat_result = statistics()
		else:
			raise ValueError("you have entered a false param. Param can only be either ALL, analysis or stat")



def analysis(tmp_filetype, tmp_path, tmp_design_file):

	if tmp_filetype == "FASTQ":
		FQanalyse = analysisFQ()
		FQanalyse.get_file(tmp_path)
		analyse_result = FQanalyse.count_read(tmp_design_file)
		print analyse_result
	elif tmp_filetype == "BAM":
		BAManalyse = analysisBAM()
		BAManalyse.get_file(self.file_path)
		analyse_result = BAManalyse.count_read(tmp_design_file)
	return analyse_result	

		



if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	config = ConfigParser.SafeConfigParser()

	parser.add_argument('--designfile',required=True,help='path to the design file, the design file is the file containing the variant')
	parser.add_argument('path', help='where the sample file are stored')
	parser.add_argument('-t','--analysistype', default='ALL', help='type of analysis to be processed')
	parser.add_argument('-f','--filetype', default='FASTQ', help='file type to process')
	args = parser.parse_args()
	config.read(args.designfile)
	print config
	hotcount_stda = stand_alone(config, args.path, args.analysistype, args.filetype)
	hotcount_stda.choose_analysis()
	print hotcount_stda.analyse_result