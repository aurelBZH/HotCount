#/usr/bin/python
from Bio import SeqIO
import glob2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqUtils
import ipdb

class HotCount(object):

		def __init__(self):
			self.design_file = ""
			self.data_table = ""
			self.analyse_results = {}
			self.statistics_results = ""
		
		def analyse(self,design_file, path, filetype):                                                                            

			print(design_file.items('mutation_design'))
			design_dict={}
			for design_values in design_file.items('mutation_design'):
				design_dict[design_values[0]] = Seq(design_values[1]) 

			if filetype == "FASTQ":
				file_extentions = patterns = ['*.fastq.*', '*.fq.*', '*.fastq','*.fq','*.FASTQ']
				for extention in file_extentions : 
					
					file_list = glob2.glob(path+extention)
					print file_list
					if len(file_list) != 0 :
						break
				if len(file_list) == 0 :	 
					raise Exception("the path is incorrect or the file extension"+ 
						"is bad,in fact the program can't find the file ")


                                
				for file in file_list:
					self.analyse_results[file] = {}
					mutation_number_file_variant = []
					with open(file, "rU") as handle:
						records = list(SeqIO.parse(handle, filetype.lower()))
          
						for name, design in design_dict.iteritems():
							mutation_number_by_var_val = 0       
							print design
							mut_number = 0
							for record in records : 
								if len(SeqUtils.nt_search(str(record.seq),str(design)))>1:
									mut_number = mut_number+1
								elif len(SeqUtils.nt_search(str(record.seq),str(design.reverse_complement())))>1: 
									mut_number = mut_number+1 
							print mut_number
							mutation_number_by_var_val= mutation_number_by_var_val + mut_number
							print mutation_number_by_var_val
							#ipdb.set_trace()
							mutation_number_file_variant.append([name,mutation_number_by_var_val])
					self.analyse_results[file] = mutation_number_file_variant                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
			
				print self.analyse_results
					
			elif filetype == "BAM":
				pass
			else : 
				raise TypeError(" the file type you'll have chosen isn't processed by the program")

		def statistics(self):
			print("why")
			print("test")
	
