#/usr/bin/python
# coding: utf-8 
import pytest
import os 
from src import hotcount_standalone


def test_writetofile(tmpdir):
	count = {'TS_L251a_21_Ampliseq_2014-11-14.fastq': {'all': 550, 'wt': 282, 'delg': 24},
			 'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': {'all': 317, 'wt': 151, 'delg': 22}}
	res_file = tmpdir.mkdir("sub").join("result.txt")
	file = hotcount_standalone.to_csv(count, str(res_file))
	with open(file, 'r') as file:
		assert file.readline() == 'sample,all,wt,delg\r\n'

def test_to_dict():
	assert diff(hotcount_standalone.to_dict("resources/testfile/test2.txt"),
				{'/home/aurelien/HotCount_project/hotcount/data/asxl1/dupG_Pos/TS_L251a_21_Ampliseq_2014-11-14.fastq': {'all': '750', 'delg': '420', 'wt': '140'},
				 "/home/aurelien/HotCount_project/hotcount/data/asxl1/dupG_Pos/TS_L1095a_055_Ampliseq_AML_201-12-22.fastq":{'all':'317','delg':'151','wt':'22'}})==0


def test_count():
	hotcount_standalone.StandAlone().count(
		{'all': 'ATCGGAGGG.*GGGTGGCCC', 'wt': 'ATCGGAGGGGGGGGTGGCCC', 'delg': 'ATCGGAGGGGGGGTGGCCC'},
		'resources/testfile/asxl1/file1/','FASTQ',True,"DNA","result.txt")
	with open("result.txt", 'r') as file:
		assert file.readline() == 'sample,all,wt,delg\r\n'
		assert file.readline() == "resources/testfile/asxl1/file1/TS_L251a_21_Ampliseq_2014-11-14.fastq, 550,282,24\r\n"
def test_stat():
	#hotcount_standalone.StandAlone().stat()
    pass
def test_all():
	pass



def diff(list1, list2):
	"""	
	docstring pour diff 
	function used to compare 2 list without redundancy 
	usable for test with known input 

	 """
	c = set(list1).union(set(list2))
	d = set(list1).intersection(set(list2))
	return len(list(c - d))