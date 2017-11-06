#! /usr/bin/env python 2.7
# coding: utf-8

import logging
#import graypy
from logging.handlers import RotatingFileHandler



logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')


file_handler = RotatingFileHandler('hotcount.log', 'a', 1000000, 1)
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


steam_handler = logging.StreamHandler()
steam_handler.setLevel(logging.INFO)
logger.addHandler(steam_handler)


# graypyhandler = graypy.GELFHandler('http://10.9.206.151/', 12201, debugging_fields=True)
# print(graypyhandler)
# logger.addHandler(graypyhandler)


