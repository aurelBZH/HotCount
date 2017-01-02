#/usr/bin/python
from Bio import SeqIO
import glob2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
class HotCount(object):

		def __init__(self):
			self.design_file = ""
			self.data_table = ""
			self.analyse_results = {}
			self.statistics_results = ""
		
		def analyse(self,design_file, path, filetype):

			print(design_file.items('mutation_design'))
			print(filetype)

				
			if filetype == "FASTQ":
				file_extentions = patterns = ['*.fastq.*', '*.fq.*', '*.fastq','*.fq','*.FASTQ']
				for extention in file_extentions : 
					file_list = glob2.glob(path+extention)
					
					if len(file_list) != 0 :
						break
				if len(file_list) == 0 :	 
					raise Exception("the path is incorrect or the file extension is bad,in fact the program can't find the file ")
				 
				for file in file_list:
					with open(file, "rU") as handle:
						for record in SeqIO.parse(handle, filetype.lower()):
							for design_values in design_file.items('mutation_design'):
								Seq(design_values[1],generic_dna)

								record.seq.find(Seq(design_values[1],generic_dna))	
							# print(record)
							pass

			elif filetype == "BAM":
				pass
			else : 
				raise TypeError(" the file type you'll have chosen isn't processed by the program")

		def statistics(self):
			print("why")
			print("hihi")
	
