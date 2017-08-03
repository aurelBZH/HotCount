# !/usr/bin/env python
# coding: utf-8

from setuptools import setup, find_packages
import python_Bio_Regexp
import hotcount
import hotcount_standalone


setup(
name ="hotcount_remastered",

packages = find_packages(),
install_requires=["alabaster==0.7.10",
"appdirs==1.4.3",
"Babel==2.4.0",
"backports.shutil-get-terminal-size==1.0.0",
"biopython==1.68",
"ConfigArgParse==0.11.0",
"configparser==3.5.0",
"decorator==4.0.10",
"docutils==0.13.1",
"enum34==1.1.6",
"glob2==0.5",
"graypy==0.2.14",
"imagesize==0.7.1",
"ipython-genutils==0.1.0",
"Jinja2==2.9.6",
"logging==0.4.9.6",
"MarkupSafe==1.0",
"multiprocessing==2.6.2.1",
"numpy==1.11.3",
"packaging==16.8",
"panda==0.3.1",
"pandas==0.20.1",
"pathlib2==2.1.0",
"pexpect==4.2.1",
"pickleshare==0.7.4",
"prompt-toolkit==1.0.9",
"ptyprocess==0.5.1",
"py==1.4.32",
"Pygments==2.1.3",
"pyparsing==2.2.0",
"pysam==0.11.2.2",
"pytest==3.0.5",
"python-dateutil==2.6.0",
"pytz==2017.2",
"requests==2.14.2",
"scipy==0.19.0",
"simplegeneric==0.8.1",
"six==1.10.0",
"snowballstemmer==1.2.1",
"Sphinx==1.6.2",
"sphinxcontrib-websupport==1.0.1",
"traitlets==4.3.1",
"typing==3.6.1",
"wcwidth==0.1.7"],

author= "Aurélien Béliard, Martin Figeac",

author_email = "aurelien.beliard@univ-lille1.fr",

description =" a software to find genomic mutation with regular expression",

long_description = open('README.md').read(),

include_package_data = True,

url = "https://github.com/aurelBZH/HotCount.git",

classifiers=[
        "Programming Language :: Python",
        "Development Status :: 1 - Planning",
        "License :: OSI Approved",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2.7",
        "Topic :: genomic",
    ],
)