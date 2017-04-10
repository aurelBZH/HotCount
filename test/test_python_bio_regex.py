import pytest
from hotCount_remastered import python_Bio_Regexp 

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

@pytest.mark.xfail
def test_parenthesis(reg):
	reg.__init__()
	assert reg.regex_reverse(r"TA(AC)CT")=="TC(CA)AT","first test on parenthesis "
	
	reg.__init__()
	assert reg.regex_reverse(r"A(TG)*C") == "C(GT)*A", "test star plus parenthesis"
	reg.__init__()
	assert reg.regex_reverse(r"(AC|BW)")=="(CA|WB)", "parenthesis plus pipe"
	reg.__init__()
	assert reg.regex_reverse(r"(AC|BW(GG(AC)(CC)))") == "(CA|WB(GG(CA)(CC)))", " imbricated parenthesis"
	reg.__init__()
	assert reg.regex_reverse(r"C(TG){5,8}") == "(GT){5,8}C"


def test_othersymbol(reg):
	reg.__init__()
	assert reg.regex_reverse(r"^ATCG") == "^GCTA", " simple test "

@pytest.fixture
def reg():
	return python_Bio_Regexp.regex_seq_finder()
	
	
	