"""
### purpose
# make sure these modules can be loaded before starting pipline
# this is a list of modules across .py scripts
# also a source for functions common across .py scripts
###

### fix
# export PYTHONPATH in .sh files to import these functions
###
"""

import os
# import sys # commented imports are not used in these funcs below, but are used in other python apps in the pipeline
# import json
# import math
# import time
import pickle
# import random
# import compiler
# import numpy as np
# import pandas as pd
from os import path as op
# from collections import OrderedDict, Counter


def fs(directory):
    return sorted([op.join(directory, f) for f in os.listdir(directory)])


def pkldump(obj, f):
    with open(f, 'wb') as o:
        pickle.dump(obj, o, protocol=pickle.HIGHEST_PROTOCOL)


def pklload(path):
    pkl = pickle.load(open(path, 'rb'))
    return pkl


def get_email_info(parentdir):
    from coadaptree import pklload
    pkl = op.join(parentdir, 'email_opts.pkl')
    return pklload(pkl)


def uni(mylist):
    return list(set(mylist))


def luni(mylist):
    return len(uni(mylist))


def makedir(directory):
    if not op.exists(directory):
        os.makedirs(directory)
    return directory


def createdirs(dirs):
    for d in dirs:
        makedir(d)
