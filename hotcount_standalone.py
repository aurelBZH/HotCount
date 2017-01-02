#/usr/bin/python

import argparse
import configargparse
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
			# self.stat_result = HtC.statistics()
			print "muhahah"
			print self.analyse_result
		elif self.analysis_type == "analysis":
			self.analyse_result = HtC.analyse(self.design_file, self.file_path, self.filetype)
		# elif self.analysis_type == "stat":
		# 	self.stat_result = HtC.statistics()
		else:
			raise ValueError("you  have entered a false param. Param can only be either ALL, analysis or stat")





		



if __name__ == '__main__':
	
	parser = configargparse.getArgumentParser()
	parser.add_argument('design_file',help='the')
	parser.add_argument('path', help='where the file are stored')
	parser.add_argument('-t','--type', default='ALL')
	parser.add_argument('-f','--filetype', default='FASTQ')
	args = parser.parse_args()
	hotcount_stda = stand_alone(args.design_file, args.path, args.type, args.filetype).choose_analysis() 