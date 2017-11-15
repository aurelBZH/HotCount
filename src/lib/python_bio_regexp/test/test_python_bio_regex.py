#/usr/bin/python
# coding: utf-8 
import pytest
from src.lib.python_bio_regexp.src import python_Bio_Regexp
import re
# from datadiff.tools import assert_equal 
# from datadiff import diff
# @pytest.mark.xfail

def test_simplecase(reg):
	
	# with pytest.raises(AssertionError) as excinfo:
	# assert reg.regex_reverse(r"AW{1,11}A.*TTSS(CG){2,22}AA.(AT)")=="(TA).AA(GC){2,22}SSTT.*AW{1,11}A"
	reg.__init__()
	assert reg.regex_reverse(r"ATCG") == "GCTA", " simple test "
	reg.__init__()
	assert reg.regex_reverse(r"A[TG]*C") == "C[TG]*A", "test star plus bracket"

	reg.__init__()
	assert reg.regex_reverse(r'AT.') == ".TA", "test point"
	reg.__init__()
	assert reg.regex_reverse(r"AT*") == "T*A", "test star" 
	reg.__init__()
	assert reg.regex_reverse(r"C.*T")=="T.*C", "test point star"
	reg.__init__()
	assert reg.regex_reverse(r"C.*T")=="T.*C", "test point star"
	# reg.__init__()
	# assert reg.regex_reverse(r"A[TG*]C") == "C [TG*]A", "test star in bracket"
	reg.__init__()
	assert reg.regex_reverse(r"AC[ATC]CC")=="CC[ATC]CA", "test bracket"

def test_bracecase(reg):
	reg.__init__()
	assert reg.regex_reverse(r"AW{1,11}A") == "AW{1,11}A", "simple brace"
	reg.__init__()
	assert reg.regex_reverse(r"AC{4,20}G{11,22}") == "G{11,22}C{4,20}A", "2 brace in one "
	reg.__init__()
	assert reg.regex_reverse(r"A[CG]{1,8}") == "[CG]{1,8}A", "brace after bracket"
	reg.__init__()

# @pytest.mark.xfail
def test_parenthesis(reg):
	reg.__init__()
	assert reg.regex_reverse(r"TA(AC)CT")=="TC(CA)AT","first test on parenthesis "
	
	reg.__init__()
	# assert reg.regex_reverse(r"A(TG)*C") == "C(GT)*A", "test star plus parenthesis"
	reg.__init__()
	assert reg.regex_reverse(r"(AC|BW)")=="(WB|CA)", "parenthesis plus pipe"
	reg.__init__()
	assert reg.regex_reverse(r"(AC|BW(GG(AC)(CC)))") == "(((CC)(CA)GG)WB|CA)", "(CA|WB(GG(CA)(CC)))imbricated parenthesis"
	reg.__init__()
	assert reg.regex_reverse(r"C(TG){5,8}") == "(GT){5,8}C"


def test_othersymbol(reg):
	reg.__init__()
	assert reg.regex_reverse(r"^ATCG") == "GCTA^", " simple test "


def test_verify_regex_simple_DNA(reg):
	with pytest.raises(Exception) as excinfo:
		reg.verify_regex("djgfkGDOLgeuz")
	excinfo.match(r"false symbol in the regular expression")
def test_verify_regex_DNA(reg):
	reg.verify_regex("ATC")
def test_verify_regex_DNA(reg):
	with pytest.raises(Exception) as excinfo:
		reg.verify_regex("ATC{AA}")
	excinfo.match(r"malformed regular expression")

def test_verify_regex_DNA(reg):
	# with pytest.raises(Exception) as excinfo:
	reg.verify_regex("ATC{1,11}")
	# excinfo.match(r"malformed regular expression")

def test_verify_regex_DNA_min(reg):
	with pytest.raises(Exception) as excinfo:
		reg.verify_regex("atc")
	excinfo.match(r"false symbol in the regular expression")

def test_verify_regex_DNA_brace(reg):
	with pytest.raises(Exception) as excinfo:
		reg.verify_regex("ATC{A,2}")
	excinfo.match(r"malformed regular expression")

def test_verify_regex_RNA(reg):
	with pytest.raises(Exception) as excinfo:
		reg.verify_regex("ATG", nuctype = "RNA")
	excinfo.match(r"false symbol in the regular expression")

def test_verify_regex_RNA(reg):
	with pytest.raises(Exception) as excinfo:
		reg.verify_regex(r"AT(ACC)CCT", nuctype = "RNA")
	excinfo.match(r"false symbol in the regular expression")	

def test_verify_regexDNA_bracket(reg):
		reg.verify_regex(r"AT[AC]CCT", nuctype = "DNA")

def test_verify_regex_RNA(reg):
	reg.verify_regex("AUG", nuctype = "RNA")

def test_find_subseq3(reg):
	with pytest.raises(AssertionError) as excinfo:
		assert reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", r"ATCT{1,12}", False, True, False)== ["ATCT{1,12}",0]


def test_regex_complement(reg):
	assert reg.regex_complement("ATCG") == "TAGC", "simple DNA complement"

def test_regex_complement_brace(reg):
	assert reg.regex_complement(r"A[TG]*C") == "T[AC]*G", "simple brace case "

def test_regex_complement_bracket(reg):
	assert reg.regex_complement(r"AW{1,11}A") == "T[AT]{1,11}T", "bracket case"

def test_regex_complement_bracket(reg):
	assert reg.regex_complement(r"AC{4,20}G{11,22}") == "TG{4,20}C{11,22}", "double bracket case"

def test_regex_complement_bracket(reg):
	assert reg.regex_complement(r"(AC|BW(GG(AC)(CC)))") == "(TG|[ACG][AT](CC(TG)(GG)))", "double parenthesis case"

def test_regex_complement_RNA(reg):	
	assert reg.regex_complement(r"AUCG", nuctype="RNA") == "UAGC", "simple RNA test"

def test_regex_reverse_complement1(reg):
	assert reg.regex_reverse_complement(r"ATCCT") == "AGGAT" 

def test_regex_reverse_complement2(reg):
	assert reg.regex_reverse_complement(r"ATC{1,10}CT") == "AGG{1,10}AT" 

def test_regex_reverse_complement3(reg):
	assert reg.regex_reverse_complement(r"AT[GA]CCT") == "AGG[CT]AT" 

def test_regex_reverse_complement4(reg):
	assert reg.regex_reverse_complement(r"AT(ACC)CCT") == "AGG(GGT)AT" 

def test_find_subseq(reg):
	pattern = reg.use_iupac(r"AW{1,10}(CG){1,10}", "DNA")

	comp_pattern = re.compile(pattern)
	assert diff(reg.find_subseq("ATCTTTTTATTTCGCGCGGGGAAA",comp_pattern,True,True, True, False, False), ['A[AT]{1,10}(CG){1,10}', 1]) == 0

def test_find_subseq2(reg):
	pattern = re.compile(r"ATC")
	assert diff(reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", pattern,True,False, True, False, False), ["ATC", 3]) == 0

def test_find_subseq2b_position_of_match(reg):
	assert diff(reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", r"ATC",False,True, False, True, False), ['ATC', 0, 8, 18]) == 0

def test_find_subseq2t(reg):
	assert diff (reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", r"ATC",False,True, False, False, True), ["ATC", True]) == 0


def test_find_subseq3(reg):
	pattern = re.compile(r"ATCT{1,12}")
	assert diff(reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", pattern,True,True, True, False, False), ["ATCT{1,12}", 2]) == 0


def test_find_subseq3b(reg):
	pattern = re.compile(r"ATCT{1,12}")
	tmp = reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", pattern,True,True, False, False, True)
	assert diff(tmp, ["ATCT{1,12}", True])== 0 	

def test_find_subseq3t(reg):
	pattern = re.compile(r"AT(CT){1,12}")
	assert diff(reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", pattern,True, True, False, False, True), ["AT(CT){1,12}", True]) == 0

def test_find_subseq4(reg):
	assert diff (reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", r"AT(CT){1,12}",False, False, True, False, False), ["AT(CT){1,12}", 2]) == 0

def test_find_subseq4b(reg):
	pattern = reg.use_iupac(r"AT(CT){1,12}", "DNA")
	comp_pattern = re.compile(pattern)
	assert diff(reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", comp_pattern,True, False, False, True, False), ["AT(CT){1,12}", 0, 8]) == 0

def test_find_subseq4t(reg):
	pattern = re.compile(r"AT(CT){1,12}")
	assert diff(reg.find_subseq("ATCTTTTTATCTCGCGCGATCGAAA", pattern,True, False, False, False, True),["AT(CT){1,12}", True]) == 0

def test_use_IUPAC_for_w(reg):
	assert reg.use_iupac(r"AW{1,11}A","DNA") == "A[AT]{1,11}A"


def test_use_IUPAC_for_w(reg):
	assert reg.use_iupac(r"RY{1,11}S","DNA") == "[AG][CT]{1,11}[CG]"

def test_create(reg):
	pattern=reg.create_pattern(r"RY{1,11}S","DNA",True )
	assert pattern.pattern == "[AG][CT]{1,11}[CG]"

# fonction utilisé pour comparer 2 tableau sans redondance 
#utilisable dans la mesure ou les test ont une entrée connu
#
def diff(list1, list2):
	"""	docstring pour diff  """
	c = set(list1).union(set(list2))
	d = set(list1).intersection(set(list2))
	return len(list(c - d))



@pytest.fixture
def reg():
	return python_Bio_Regexp.regex_seq_finder()



