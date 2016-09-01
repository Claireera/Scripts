# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 15:13:59 2016

@author: claire
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 14:09:51 2016

@author: claire

Plot a box plot of pEak to Peak a s a function of Azimuth
"""

import h5py
import numpy as np
from Events_Selections import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *
from Instance_EQPeaktoPeak import *
from Signal_Wave_Picking import *
import cPickle
import copy
import pandas as pd

MlMin= 5
MlMax = 6
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'

EQ = 1
Station = '1'
#Read File with EQ results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f56 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V2_Ml_%s_%s_distmax_%s.hdf5'%(str(MlMin),str(MlMax),str(DistMax)),'r+')
g= f56.keys()
h = copy.copy(g)

DPeaktoPeakMean ={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DPGA = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DPGV = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DArias ={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DML={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DAz = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
Drdist = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}

for i in h: 
    k=f56[i].keys()
    for Station in k: 
        for Component in ['Z','R','T','E','N','H']:
            if f56[i][Station]['valid_{}'.format(Component)][0]==True and f56[i][Station]['Satured_{}'.format(Component)][0]==False : 
                #print Component, f[i][Station]['PtP_{}'.format(Component)][:],'sature', f[i][Station]['SNR_{}'.format(Component)][:]
                
                DML[Station][Component].append(f56[i].attrs['Ml'])
                DAz[Station][Component].append(f56[i][Station].attrs['Az'])
                Drdist[Station][Component].append(f56[i][Station].attrs['RDist'])
                DPGA[Station][Component].append(f56[i][Station]['PGA_{}'.format(Component)][0])
                DPGV[Station][Component].append(f56[i][Station]['PGV_{}'.format(Component)][0])
                DArias[Station][Component].append(f56[i][Station]['PGV_{}'.format(Component)][0])
                DPeaktoPeakMean[Station][Component].append(np.mean(np.array(f56[i][Station]['PtP_{}'.format(Component)][0:3])))
                
                
                
#Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                
LStations = ['1','2','3','4','5','6','7']
#component 
LComponents = ['N','R','E','T','Z','H']

for st in LStations : 
    Station = 'St'+st
    fig,  axes = plt.subplots(2, 3, subplot_kw=dict(projection='polar'))
    #fig.subplots_adjust(hspace=0.15, wspace=0.3)
  
    a =0
    for Component in LComponents :
        
        #subplot number
        if a<=2 :
            ax=axes[0,a]
        else :
            ax=axes[1,a-3]
        
        a +=1
        #form tuple angle, Lparam
        LParam = DPeaktoPeakMean[Station][Component]
        Langle = DAz[Station][Component]
        Ltuple =[]
        for i in xrange(len(Langle)):
            Ltuple.append((Langle[i],LParam[i]))
        
        angle = math.radians(10.)
        patches = math.radians(360.)/angle
        theta = np.arange(0,math.radians(360.),angle)
        count = [[] for x in xrange(int(patches))]

        countbis = [0.]*int(patches)
        for angi, Param in Ltuple :
            print 'param', Param, angi
            if angi<0:
                #angle is given between -180;180 
                ang=360+angi
            else:
                ang=angi
                
            item=math.radians(ang)
            temp = int((item - item%angle)/angle)
            count[temp].append(Param)