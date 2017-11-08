#/usr/bin/python
# coding: utf-8 
import pytest
from src import python_Bio_Regexp
from src import hotcount
from src import hotcount_standalone

# def test_get_file():
# 	cls = hotcount.AnalysisFQ()
# 	with pytest.raises(Exception) as excinfo:
# 		print (type(excinfo))
# 		cls.get_file("/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1")
# 	assert excinfo.value == ""

	#check if something raise an error
# def test_get_file2():
# 	cls=hotcount.AnalysisFQ()
# 	tab = cls.get_file("/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1/file/")
# 	print (tab)
# 	assert tab == ['/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1/file/TS_L1095a_055_Ampliseq_AML_201-12-22.fastq', '/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1/file/TS_L1093a_57_Ampliseq_2015-01-13_2.fastq', '/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1/file/TS_L1106a_62_Ampliseq_2015-01-13.fastq', '/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1/file/TS_L1103a_60_Ampliseq_2015-01-13_2.fastq', '/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1/file/TS_L1101a_58_Ampliseq_2015-01-13_2.fastq', '/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1/file/TS_L1102a_59_Ampliseq_2015-01-13_2.fastq', '/home/aurelien/HotCount_project/hotCount_remastered/resources/testfile/asxl1/file/TS_L251a_21_Ampliseq_2014-11-14.fastq']
#
def test_count():
	FQanalyse = hotcount.AnalysisFQ()
	FQanalyse.get_file("resources/testfile/asxl1/file1/")
	analyse_result = FQanalyse.count_read({'all': 'ATCGGAGGG.*GGGTGGCCC', 'wt': 'ATCGGAGGGGGGGGTGGCCC', 'delg': 'ATCGGAGGGGGGGTGGCCC'},"DNA",True)
	assert diff(analyse_result, {'resources/testfile/asxl1/file1/TS_L251a_21_Ampliseq_2014-11-14.fastq': {'all': 550, 'delg': 24, 'wt': 282}}) == 0

def test_count2():
	FQanalyse = hotcount.AnalysisFQ()
	FQanalyse.get_file("resources/testfile/asxl1/file1gz/")
	analyse_result = FQanalyse.count_read({'all': 'ATCGGAGGG.*GGGTGGCCC', 'wt': 'ATCGGAGGGGGGGGTGGCCC', 'delg': 'ATCGGAGGGGGGGTGGCCC'},"DNA",False)
	assert diff(analyse_result, {'resources/testfile/asxl1/file1gz/TS_L251a_21_Ampliseq_2014-11-14.fastq.gz': {'all': 550, 'delg': 24, 'wt': 282}}) == 0

def test_all():
	pass


def test_stat(stat):
	res = diff(stat.create_contingency_table(),{'wt': {'TS_L251a_21_Ampliseq_2014-11-14.fastq': [282, 550], 'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': [151, 317]}, 'delg': {'TS_L251a_21_Ampliseq_2014-11-14.fastq': [240, 550], 'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': [22, 317]}}
)
	assert res == 0

def test_stat2(stat):
	stat.create_contingency_table()
	assert diff(stat.apply_fisher_test(), {'wt': {'TS_L251a_21_Ampliseq_2014-11-14.fastq': {'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': 0.29613043222614138}, 'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': {'TS_L251a_21_Ampliseq_2014-11-14.fastq': 0.74487252873289977}}, 'delg': {'TS_L251a_21_Ampliseq_2014-11-14.fastq': {'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': 2.9284053112022063e-21}, 'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': {'TS_L251a_21_Ampliseq_2014-11-14.fastq': 1.0}}}
) == 0


def test_stat3(stat):
 	stat.create_contingency_table()
 	stat.apply_fisher_test()
	assert diff(stat.check_positiv_sample(), {'wt': [['TS_L251a_21_Ampliseq_2014-11-14.fastq', 1, 0.29613043222614138], ['TS_L1095a_055_Ampliseq_AML_201-12-22.fastq', 1, 0.74487252873289977]], 'delg': [['TS_L251a_21_Ampliseq_2014-11-14.fastq', 0, 2.9284053112022063e-21], ['TS_L1095a_055_Ampliseq_AML_201-12-22.fastq', 1, 1.0]]}) == 0

#ne test pas tout ce qu'il devrait !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def test_all_stat():
 	count = {'TS_L251a_21_Ampliseq_2014-11-14.fastq': {'all': 550, 'wt': 282, 'delg': 24}, 'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': {'all': 317, 'wt': 151, 'delg': 22}, 'TS_L1095a_055_Ampliseq_AML_202-11-23.fastq': {'all': 317, 'wt': 180, 'delg': 130},  'TS_L1095a_055_Ampliseq_AML_203-8-28.fastq': {'all': 417, 'wt': 151, 'delg': 22}}
 	assert 	diff(hotcount_standalone.global_stat(count, 0.01, 6,"all",100, ['delg', 'wt']),{'wt': [['TS_L251a_21_Ampliseq_2014-11-14.fastq', [0.29613043222614138]], ['TS_L1095a_055_Ampliseq_AML_201-12-22.fastq', [0.74487252873289977]]], 'delg': [['TS_L251a_21_Ampliseq_2014-11-14.fastq',[2.9284053112022063e-21]], ['TS_L1095a_055_Ampliseq_AML_201-12-22.fastq',[1.0]]]}) == 0


# function used to compare 2 list without redundancy 
#usable for test with known input 

def diff(list1, list2):
	"""	
	docstring pour diff 
	function used to compare 2 list without redundancy 
	usable for test with known input 

	 """
	c = set(list1).union(set(list2))
	d = set(list1).intersection(set(list2))
	return len(list(c - d))



@pytest.fixture
def stat():
	count = {'TS_L251a_21_Ampliseq_2014-11-14.fastq': {'all': 550, 'wt': 282, 'delg': 240}, 'TS_L1095a_055_Ampliseq_AML_201-12-22.fastq': {'all': 317, 'wt': 151, 'delg': 22}}

	return hotcount.statistics(count,0.01, 6, "all",100, ["wt","delg"])
