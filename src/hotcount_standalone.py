#! /usr/bin/env python 2.7
# coding: utf-8
from __future__ import print_function

import ConfigParser
import argparse
import csv
from hotcount import *
import gzip

version = 1.0


class StandAlone(object):
    """

    """
    def __init__(self):
        self.analyse_result = ""
        self.stat_result = ""

	

    def count(self, design_dict, path, filetype, IUPAC, nuctype,output_file="default"):
        """
        a function making the count treatment only
        :param design_dict: file containing the design of the  diferent regex to count
        :param path: path to the directory where the file are stored
        :param filetype :file type
        :param output_file : number of positiv sample chosen by the user
        :type design_dict: dictionary
        :type path : string
        :type filetype: string
        :type output_file: string
        """
        if output_file!=None:
            self.analyse_result = analysis(filetype, path, design_dict,nuctype, IUPAC)
            try:
                op_file = open(output_file, "w")
            except Exception as e:
                raise Exception("there is a problem with the file opening")
            to_csv(self.analyse_result, output_file)
        else:
            self.analyse_result = analysis(filetype, path, design_dict,nuctype, IUPAC)
#            csv_analysis = to_csv(self.analyse_result, output_file)
            logger.info(self.analyse_result)


    def ALL(self,design_dict, path, filetype, pvalue, mutation, sample, controle, depth, IUPAC,nuctype, output_file="default" ):
        """
        a function to make count and statistic treatment
        and display the result .
        :param design_dict: file containing the design of the  diferent regex to count
        :param path: path to the directory where the file are stored
        :param filetype : number of positiv sample chosen by the user
        :param pvalue: pvalue to choose for analysis
        :param mutation: list of mutation to study
        :param sample : number of positiv sample chosen by the user
        :param controle: controle value for statistics
        :param output_file : number of positiv sample chosen by the user
        :type design_dict: dictionary
        :type path : string
        :type filetype : string
        :type pvalue: float
        :type mutation : string
        :type sample : float
        :type controle : string
        :type output_file: string
        """
        if output_file!=None:
            try:
                op_file = open(output_file, "w")
            except Exception as e:
                raise Exception("there is a problem with the file opening")
            self.analyse_result = analysis(filetype, path, design_dict, nuctype, IUPAC)
            to_csv(self.analyse_result, output_file)
            self.stat_result = global_stat(self.analyse_result, pvalue, sample, controle, depth, mutation)
            for k in xrange(0,int(sample)):
                for i,j in self.stat_result.iteritems():
                    if j==k:
                        logger.debug(str(i)+','+str(j))
                        op_file.write(str(i)+','+str(j))
            op_file.close()
        else:
            self.analyse_result = analysis(filetype, path, design_dict,nuctype, IUPAC)

            self.stat_result = global_stat(self.analyse_result, pvalue, sample, controle, depth, mutation)
            for k in xrange(0,int(sample)):
                for i,j in self.stat_result.iteritems():
                    if j==k:
                        logger.info (i,j)

    def stat(self, countfile, pvalue, sample, mutation, controle,depth, output_file="default"):
        """
        a function to make count and statistic treatment
        and display the result
        :param countfile: file containing the result of the  diferent count
        :param pvalue: pvalue to choose for analysis
        :param sample : number of positiv sample chosen by the user
        :param mutation: list of mutation to study
        :param controle: controle value for statistics
        :param output_file : number of positiv sample chosen by the user
        :type countfile: dictionary
        :type pvalue: float
        :type sample : float
        :type mutation : string
        :type controle : string
        :type output_file: string
        """
        file_ext =[".*\.txt",".*\.csv"]
        for i in file_ext:
            if re.match(i, countfile):
                countvalue = to_dict(countfile)
        if output_file!=None:
            try:
                op_file = open(output_file, "w")
            except Exception as e:
                raise Exception("there is a problem with the file opening")
            self.stat_result =global_stat(countvalue, pvalue, sample, controle, depth, mutation)
            x=0
            logger.debug(sample)
            for mut,result in self.stat_result.items():
                op_file.write("results for the "+str(mut)+" mutation vs all")
                logger.info("results for the "+str(mut)+" mutation vs all")
                positif_number=int(len(result[1][1])) if int(len(result[1][1]))<=int(sample) else int(sample)
                for i in range(0, int(len(result[1][1]))):
                    op_file.write("Who is significant with only  " + str(
                        i+1) + " positive  sample  in this library (p<=" + str(pvalue) + ")")
                    logger.info("Who is significant with only " + str(
                        i+1) + " positive  sample  in this library (p<=" + str(pvalue) + ")")
                    a=False
                    count_significativ_sample=0
                    for result_by_file in result:
                        if result_by_file[1][i]<=pvalue:
                            count_significativ_sample+=1
                            a=True
                            op_file.write(
                               "the sample" + str(result_by_file[0]) + " is statisticaly significant with a p-value of "+str(result_by_file[1][i]))
                            logger.info(
                               "the sample" + str(result_by_file[0]) + " is statisticaly significant with a p-value of " + str(
                                   result_by_file[1][i]))
                    if count_significativ_sample !=i:
                        op_file.write("the significativ sample number don't match the hypothesis")
                        logger.info("the significativ sample number don't match the hypothesis")
                    if(a==False):
                       op_file.write("none")
                       logger.info("none")



def analysis(tmp_filetype, tmp_path,tmp_design_dict,nuctype, IUPAC):
    """
    a function regrouping the count analysis based on python bio regex

    :param tmp_filetype: input file type
    :param tmp_path: path of input files
    :param tmp_design_dict: design_dict containing design
    :type tmp_filetype: string
    :type tmp_path: string
    :type tmp_design_dict: string
    :return analyse_result: result of the count analysis
    :rtype analyse_result: dictionary

    """
    if tmp_filetype == "FASTQ":
        FQanalyse = AnalysisFQ()
        FQanalyse.get_file(tmp_path)
        analyse_result = FQanalyse.count_read(tmp_design_dict,nuctype, IUPAC)
    elif tmp_filetype == "BAM":
        BAManalyse = AnalysisBAM()
        BAManalyse.get_file(tmp_path)
        analyse_result = BAManalyse.count_read(tmp_design_dict,nuctype, IUPAC)
    return analyse_result


def global_stat(analyse_result, pvalue, sample, controle,depth, mutation):
    """
    a function to regroupe statistic treatment
    :param analyse_result: result of the mutation count
    :param pvalue: pvalue chosebn by the user by default 0.01
    :param sample : number of positiv sample chosen by the user
    :param controle: name of the controle string
    :param mutation: name of the mutation
    :type analyse_result: dictionary
    :type pvalue: float
    :type sample : string
    :type controle: string
    :type mutation: stringr
    :return results: results of the statistic analysis
    :rtype results: dictionary
    """
    statistic = statistics(analyse_result, pvalue, sample, controle, depth, mutation)
    statistic.create_contingency_table()
    fisher = statistic.apply_fisher_test()
    results = statistic.check_positiv_sample()
    return results


def to_csv(dicti, output_file ="default"):
    """
    a simple function to write results in csv format

    :param dicti: a dictionary of statistics result
    :param output_file: the output file to write the stat result
    :type dicti: dictionary
    :type output_file: string
    """
    cpt=0
    #if output_file == "default" or output_file == None:
     #   output_file = "result.txt"
    with open(output_file, 'wb') as f:
        for i,j in dicti.iteritems():
            w = csv.DictWriter(f,j)
            if cpt == 0 :
                f.write("sample,")
                w.writeheader()
            f.write(i+', ')
            w.writerow(j)
            cpt+=1
    return output_file


def to_dict(file):
    """
    a simple function to return a dictionnary from a csv file
    :param file: file with count result in csv format
    :type string: file name in string format
    :return result_dict: values from the treatment
    :rtype result_dict: dictionary
    """
    result_dict={}
    with open(file, 'rb') as f:
        try:
            reader = csv.DictReader(f)
            for row in reader:
                val = row["sample"]
                result_dict[val] = {}
                for key,value in row.iteritems():
                    if key != "sample":
                        result_dict[val][key]=int(value)
        except Exception as e:
            raise e("there is a problem with input file ")

    return result_dict


def show_mutation(config):
    """
    log function to show the mutation the software is working on
    :param config: an object config parser containing the diferent mutation
    :type config: configparser object
    """
    design_dict = dict(config.items('mutation_design'))
    logger.info("show mutation")
    for name, design in design_dict.iteritems():
        logger.info("%(name_v)s,%(design_v)s"%{"name_v":name,"design_v":design})

def show_arg(args):
    for i in args.iteritems():
        logger.info(i)

def init():
    """
    function for initialisation and parsing of
    parrameters

    """

    config = ConfigParser.SafeConfigParser()
    parser = argparse.ArgumentParser()

    sub = parser.add_subparsers(help="type of analysis",dest="analysis_type")
    all_parser = sub.add_parser("all")
    all_parser.add_argument('--designfile',required=True,help='path to the design file, the design file is the file containing the variant')
    all_parser.add_argument('--path',required=True, help='where the sample file are stored')
    all_parser.add_argument('-f','--filetype', default='FASTQ', help='file type to process')
    all_parser.add_argument("-o","--output", help="result file")
    all_parser.add_argument("-p","--pvalue", default=0.001, help="pvalue", type=int)
    all_parser.add_argument("-s", "--sample", default=6, help="max number of positive samples")
    all_parser.add_argument("-m", "--mutation", required= True, nargs='*', help=" mutation to be analysed in a string, separated by a ',' " )
    all_parser.add_argument("-a","--controle",default="all", help="controle regex")
    all_parser.add_argument("-d", "--depth", default=100, help="controle regex",type=int)
    all_parser.add_argument("-I", "--IUPAC", help="use IUPAC data",action='store_true')
    all_parser.add_argument("-n","--nuctype",default="DNA", help="type of nucleotide")

    count_parser = sub.add_parser("count")
    count_parser.add_argument('--designfile',required=True,help='path to the design file, the design file is the file containing the variant')
    count_parser.add_argument('--path',required=True, help='where the sample file are stored')
    count_parser.add_argument('-f','--filetype', default='FASTQ', help='file type to process')
    count_parser.add_argument("-o","--output", help="result file")
    count_parser.add_argument("-I", "--IUPAC", help="use IUPAC data",action='store_true')
    count_parser.add_argument("-n","--nuctype",default="DNA", help="type of nucleotide")


    stat_parser = sub.add_parser("stat")
    stat_parser.add_argument("-c","--countfile", required=True, help="file with count result in csv")
    stat_parser.add_argument("-p","--pvalue", default=0.001, help="pvalue", type=int)
    stat_parser.add_argument("-s", "--sample", default=6, help="max number of positive samples")
    stat_parser.add_argument("-m", "--mutation", required= True,nargs='+', help=" mutation to be analysed in a string, separated by a ',' " )
    stat_parser.add_argument("-o","--output", help="result file")
    stat_parser.add_argument("-a","--controle", default="all", help="controle regex")
    stat_parser.add_argument("-d", "--depth", default=100, help="controle regex",type=int)


    args = vars(parser.parse_args())

    if "designfile" in args:
        config.read(args["designfile"])
        show_mutation(config)
        design_dict=dict(config.items('mutation_design'))

    show_arg(args)
    hotcount_stda = StandAlone()
    analysis_type = args["analysis_type"]
    logger.info("analysis begin in %s mode" %analysis_type)
    if args["analysis_type"]=="all":
        hotcount_stda.ALL(design_dict,args["path"],args["filetype"], args["pvalue"], args["mutation"], args["sample"], args["controle"],args["depth"], args["IUPAC"],args["nuctype"], args["output"])
    elif args["analysis_type"]=="count":
        hotcount_stda.count(design_dict, args["path"], args["filetype"], args["IUPAC"],args["nuctype"], args["output"])
    else :
        hotcount_stda.stat(args["countfile"], args["pvalue"], args["sample"], args["mutation"], args["controle"], args["depth"], args["output"])


if __name__ == '__main__':
    init()