# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 15:14:12 2016

@author: claire
"""

import h5py
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import matplotlib.cm  as cm 
import matplotlib as mpl
MlMin= 3
MlMax = 3.1
DistMin = 0
DistMax = 50
jJultoplot = 206
f = h5py.File('/home/claire/PHD/Working/Results/Spectrum/Parallel_Spectrum_Ml_3_3.1_distmax_50.hdf5','r+')

g= f.keys()
i =0
for EQ in g :
    for St in ['1','2','3','4','5','6','7'] :
        if f[EQ]['St{}'.format(St)]['valid_N'][0]==False :
            print EQ,f[EQ]['St{}'.format(St)]['SNR_N'][0]
            i +=1
        