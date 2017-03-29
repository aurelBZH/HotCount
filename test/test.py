import pytest
from python_bio_regex import regex_seq_finder 

def test_python_bio_regex():
	reg = regex_seq_finder(object)
	assert reg.regex_complement(r"AW{1,11}A.*TTSS(CG){2,22}AA.(AT)")=="(AT).AA(CG){2,22}SSTT.*AW{1,11}A"