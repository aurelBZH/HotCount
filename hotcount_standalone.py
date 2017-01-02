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
		HtC = HotCount()
		if self.analysis_type == "ALL":
			self.analyse_result =HtC.analyse(self.design_file, self.file_path, self.filetype)
			self.stat_result = HtC.statistics()
			print self.stat_result
			print "muhahah"
		elif self.analysis_type == "analysis":
			self.analyse_result = HtC.analyse(self.design_file, self.file_path, self.filetype)
		elif self.analysis_type == "stat":
			self.stat_result = HtC.statistics()
		else:
			raise ValueError("you have entered a false param. Param can only be either ALL, analysis or stat")





		



if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	config = ConfigParser.SafeConfigParser()

	parser.add_argument('--designfile',required=True,help='path to the design file, the design file is the file containing the variant')
	parser.add_argument('path', help='where the sample file are stored')
	parser.add_argument('-t','--type', default='ALL', help='type of analysis to be processed')
	parser.add_argument('-f','--filetype', default='FASTQ', help='file type to process')
	args = parser.parse_args()
	config.read(args.designfile)
	hotcount_stda = stand_alone(config, args.path, args.type, args.filetype).choose_analysis() 