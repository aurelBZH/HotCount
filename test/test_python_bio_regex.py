import pytest
from hotCount_remastered import python_Bio_Regexp 

def test_python_bio_regex():
	reg = python_Bio_Regexp.regex_seq_finder()
	assert reg.regex_reverse(r"AW{1,11}A.*TTSS(CG){2,22}AA.(AT)")=="(AT).AA(CG){2,22}SSTT.*AW{1,11}A"
	reg.__init__()
	assert reg.regex_reverse(r"ATCG")=="GCTA"
	reg.__init__()
	assert reg.regex_reverse(r"[ATC]")=="[ATC]"
	# reg.__init__()
	# assert reg.regex_reverse(r"(AT)")=="(TA)"
	reg.__init__()
	assert reg.regex_reverse(r"A.*T")=="T.*A"
	reg.__init__()
	assert reg.regex_reverse(r"AW{1,11}A")=="AW{1,11}A"
	