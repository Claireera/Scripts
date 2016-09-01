# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 08:51:32 2016

@author: claire

Crop list to get samples of events
"""
import json
import os
import cPickle


MlMin= 3
MlMax = 4
DistMin = 0
DistMax = 200

with open('/home/claire/PHD/Working/Data/List_EQ_St_Comp_Freq_Ml_%s_%s_distmax_%s.txt'%(MlMin,MlMax,DistMax)) as f:
    LEqStCompFreq =json.load(f)

Listcrop =  LEqStCompFreq[17886:17952][:]
Outputfile = open('/home/claire/PHD/Working/Data/List_EQ_St_Comp_Freq_Ml_%s_%s_distmax_%s_Test.txt'%(MlMin,MlMax,DistMax), mode = 'w+')
json.dump(Listcrop,Outputfile)
Outputfile.close()
f.close()