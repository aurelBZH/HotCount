#/usr/bin/python
import argparse
import pysam

def research_mut_by_pos(poolA, poolB, output_file):

	poolA_object = pysam.AlignmentFile(poolA, "rb")
	pysam.index(poolA)
	poolB_object = pysam.AlignmentFile(poolB, "rb")
	pysam.index(poolB)
	poolA_mut_list = {"chr1":["115258747","G","A"], "chr2":["25457242","G","A"], "chr4":["106190801","T","C"],"chr7":["148512040","T","G"],"chr4":["106182956","CATTTGCAAAACCTGTCCACTCT",""],"chr5":["170837547","TCTG"]}
	poolB_mut_list = {"chr20":["31022784","C","T"], "chr12":["11905495", "C","T"], "chr2":["209113113","C","T"], "chr15":["90631934","G","A"], "chr21":["36164565","AC",""], "chr21":["36259198", "TCCT",""],"chr4":["106155861","TCAG",""], "chr4":["106190855","G","T"], "chr17":["7578406","G","A"] }
	mut_dict={}
	of = open(output_file, "w")
	for key, value in poolA_mut_list.iteritems():
		for read in poolA_object.fetch(key, int(value[0])-30, int(value[0])+30):
			if key+value[0] in mut_dict:
				mut_dict[key+value[0]].append([read.cigar, read.seq])
			else:
				mut_dict[key+value[0]] = [read.cigar, read.seq]


	of.write(str(mut_dict))
	# for key value in mut_d ict.iteritems():
# 
		
	of.close()		


	# print poolA
	# print poolB


	poolA_object.close()		

	poolB_object.close()




if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument("--poolA")
	parser.add_argument("--poolB")
	parser.add_argument("--output_file")
	args = parser.parse_args()
	research_mut_by_pos(args.poolA, args.poolB, args.output_file)