# System tools
import sys
import os
import time
import subprocess
from itertools import combinations as combi
from collections import defaultdict
import glob
from datetime import datetime
import re
import json
import random
print("system tools imported")

# Basic tools
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
print("basic tools imported")

# Chemistry tools
import ase
from ase.io import read
from molmass import Formula
from ase import io
#import dscribe
#from dscribe.descriptors import SOAP, MBTR, CoulombMatrix
print("chemistry tools imported")

# ML tools
import sklearn
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_squared_log_error, median_absolute_error
from sklearn.linear_model import LinearRegression
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C
from scipy.stats import qmc, norm
from scipy.optimize import minimize
print("ML  imported")

###################################################################
#GETTING INFORMATIONS ABOUT UFF TO USE AS INITIAL INFORMATION
_uff='./data/uff_tu'
with open(_uff,'r') as f:
    lines = f.readlines()
uff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    uff.append(words)
uff_table = pd.DataFrame(uff, columns=['atoms', 'ri', 'phi', 'xi', 'di', 'psi', 'zmm', 'vsp3', 'vsp2', 'chi', 'nc','mass'])
#grouping systems needed for the separation of atoms
#checking whether the atom is a group_6 or not
g6 = ['O','S','Se','Te','Po']


_angleff='data/angle.ff'
with open(_angleff,'r') as f:
    lines = f.readlines()
angle_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    angle_ff.append(words)
    
_torsionff='data/torsion.ff'
with open(_torsionff,'r') as f:
    lines = f.readlines()
torsion_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    torsion_ff.append(words)
    
_inversionff='data/inversion.ff'
with open(_inversionff,'r') as f:
    lines = f.readlines()
inversion_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    inversion_ff.append(words)

print('all modules and data extracted')
