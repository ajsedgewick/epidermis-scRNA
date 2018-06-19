from __future__ import  division
import numpy as np
import pandas as pd
import os

#### GENERAL PURPOSE PARSERS  ######################
def readListFromFile(fi):
    """
    Load the lines of file fi as a list
    """
    f = open(fi , 'r' )
    lines_raw = f.readlines()
    f.close()
    lines = [ elem.strip('\n')  for elem in lines_raw ]
    return lines
