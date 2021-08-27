#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 20:56:02 2018

@author: davidmurphy
"""


import numpy as np
from sys import argv
import data_processing.precalc_tools as prc
import data_processing.data_tools as dtl
from classes.runstruct import ChromStruct, root_dir

#%%
cst = ChromStruct('chr1', bdir='std_split')

