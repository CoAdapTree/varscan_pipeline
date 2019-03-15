### purpose
# make sure these modules can be loaded before starting pipline
# this is a list of modules across .py scripts
# also a source for functions common across .py scripts
###

### fix
# export PYTHONPATH in .sh files to import these functions
###

import os
import sys
import json
import math
import time
import pickle
import random
import compiler
#import numpy as np
import pandas as pd
from os import path as op
from collections import OrderedDict, Counter

def fs(DIR):
    return sorted([op.join(DIR,f) for f in os.listdir(DIR)])
def pkldump(obj,f):
    with open(f,'wb') as o:
        pickle.dump(obj,o,protocol=pickle.HIGHEST_PROTOCOL)
def pklload(path):
    pkl = pickle.load(open(path,'rb'))
    return pkl
def uni(mylist):
    return list(set(mylist))
def luni(mylist):
    return len(uni(mylist))
def makedir(DIR):
    if not op.exists(DIR):
        os.makedirs(DIR)
    return DIR
def createdirs(dirs):
    for d in dirs:
        makedir(d)
