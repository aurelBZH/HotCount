# HotCount
a programme to seach mutation on fastq or bam file without align the genome.
this programe is based on a library i built, named pythonBioRegex. This library allow to use regex in DNA RNA.you can make the reverse complement of a regex find a subseq in a seq .

installation

the best method to install our software is write this command line in the terminal:

     pip install git+https://github.com/aurelBZH/HotCount.git

you can also download the source in the github repository, open a terminal, go to the directory where the software is downloaded
and try to write :
    python setup.py install

but it can have somme bug because python easy-install is deprecated.

functionnement:

the software has 3 main function :
    -sequence count (count)
    -statistic on sequence count (stat)
    - both combined (all)


sequence count has 4 parameters:
    - --designfile file containing the regex design
    - --path for the sequence file
    - -f --filetype type of sequence file
    - -o -- output output file
 ex : python hotcount_standalone.py  count --designfile ~/HotCount_project/designfile/asxl1ini.txt --path ~/HotCount_project/hotcount/data/asxl1/dupG_Pos/ - o outputfile.txt

statistics has 6 parameters:
    - -c --countfil file with count result in csv format
    - -p --pvalue PVALUE
    - -s --sample SAMPLE max number of positive samples
    - -m --mutation MUTATION mutation to be analysed in a string, separated by a ','
    - -o OUTPUT, --output OUTPUT result file
    - -a CONTROLE, --controle CONTROLE controle regex

ex : python hotcount_standalone.py  stat  -m "mut1,mut2" -o result.txt -c ../file/to/count.txt

all has 7 parameter
    - --designfile file containing the regex design
    - --path for the sequence file
    - -f --filetype type of sequence file
    - -o -- output output file
        - -p --pvalue PVALUE
    - -s --sample SAMPLE max number of positive samples
    - -m --mutation MUTATION mutation to be analysed in a string, separated by a ','
    - -a CONTROLE, --controle CONTROLE controle regex
ex : python hotcount_standalone.py  all --designfile ~/file/where/the/design.txt --path ~/path/for/sequence/file/ -m "mut1,mut2" -o result.txt
