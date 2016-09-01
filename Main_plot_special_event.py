# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:02:01 2016

@author: claire
"""

""" Main plot special event """
from  Events_caracterisation import *
from Events_Selections import *
from Plot_results import * 
from Plot_Waveforms import *

LStations= ['1','3','5']
latevent= 22.940
longevent = 120.6
depth = 23
Year ='16'
jJul = '036'

Second = 57* 60+30
Hour= '19'
#Plot_Merging_2hours_Seismograms(Year, jJul, '05', '06', LStations)
PlotallRot(Year,jJul, Second, Hour, latevent,longevent,depth,LStations)
    