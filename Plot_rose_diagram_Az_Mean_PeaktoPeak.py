# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 14:09:51 2016

@author: claire

Plot a rose diagram of the mean values of the Mean of Peak to Peak per azimuth thershold
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



"""########################################### 
            ML 5-6
##############################################""" 
            
MlMin= 5
MlMax = 6
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'
#Read File with EQ results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f56 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s.hdf5'%(str(MlMin),str(MlMax),str(DistMax)),'r+')
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
        L=[]
        for Component in ['Z','E','N']:
           if f56[i][Station]['Exist_{}'.format(Component)][0]==False or f56[i][Station]['valid_{}'.format(Component)][0]==False or f56[i][Station]['Satured_{}'.format(Component)][0]==True: 
               L.append(False)
        if False in L : 
            continue
        else:
            for Component in ['Z','R','T','E','N','H']:
           
                #print Component, f[i][Station]['PtP_{}'.format(Component)][:],'sature', f[i][Station]['SNR_{}'.format(Component)][:]
                
                DML[Station][Component].append(f56[i].attrs['Ml'])
                DAz[Station][Component].append(f56[i][Station].attrs['Az'])
                Drdist[Station][Component].append(f56[i][Station].attrs['RDist'])
                DPGA[Station][Component].append(f56[i][Station]['PGA_{}'.format(Component)][0])
                DPGV[Station][Component].append(f56[i][Station]['PGV_{}'.format(Component)][0])
                DArias[Station][Component].append(f56[i][Station]['PGV_{}'.format(Component)][0])
                DPeaktoPeakMean[Station][Component].append(np.mean(np.array(f56[i][Station]['PtP_{}'.format(Component)][0:3])))
                
"""########################################### 
            ML 6-7
##############################################"""   

MlMin= 6
MlMax = 7
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'              
#Read File with EQ results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f67 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s.hdf5'%(str(MlMin),str(MlMax),str(DistMax)),'r+')
g= f67.keys()
h = copy.copy(g)            
for i in h: 
    k=f67[i].keys()
    
    for Station in k: 
        L=[]
        for Component in ['Z','E','N']:
           if f67[i][Station]['Exist_{}'.format(Component)][0]==False or f67[i][Station]['valid_{}'.format(Component)][0]==False or f67[i][Station]['Satured_{}'.format(Component)][0]==True : 
               L.append(False)
        if False in L : 
            continue
        else:
            for Component in ['Z','R','T','E','N','H']:
    
                #print Component, f[i][Station]['PtP_{}'.format(Component)][:an],'sature', f[i][Station]['SNR_{}'.format(Component)][:]
                
                DML[Station][Component].append(f67[i].attrs['Ml'])
                DAz[Station][Component].append(f67[i][Station].attrs['Az'])
                Drdist[Station][Component].append(f67[i][Station].attrs['RDist'])
                DPGA[Station][Component].append(f67[i][Station]['PGA_{}'.format(Component)][0])
                DPGV[Station][Component].append(f67[i][Station]['PGV_{}'.format(Component)][0])
                DArias[Station][Component].append(f67[i][Station]['PGV_{}'.format(Component)][0])
                DPeaktoPeakMean[Station][Component].append(np.mean(np.array(f67[i][Station]['PtP_{}'.format(Component)][0:3])))
                
#Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                
LStations = ['1','2','3','4','5','6','7']
#component 
LComponents = ['N','R','E','T','Z','H']

for st in LStations : 
    Station = 'St'+st
    fig,  axes = plt.subplots(2, 3, subplot_kw=dict(projection='polar'))
    fig.subplots_adjust(hspace=0.1, wspace=0.3)
  
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
         
            if angi<0:
                #angle is given between -180;180 
                ang=360+angi
            else:
                ang=angi
                
            item=math.radians(ang)
            temp = int((item - item%angle)/angle)
            count[temp].append(Param)
 
        for i in xrange(len(count)):
            countbis[i]=np.nanmean(np.array(count[i]),axis =0)*1000
            #print count[i], 'counti',np.mean(np.array(count[i]),axis =0)
            
         #width of the cololum as a function of the thershold width~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        width = angle * np.ones(patches)
        
        rmax =  np.nanmax(np.array(countbis))
        print 'rmax', np.nanmax(np.array(countbis))
        rmin = min(countbis)
        #set axis label, grid ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        ax.set_rlim(rmin,rmax)
        ax.set_theta_offset(np.pi/2)
        ax.set_thetagrids(np.arange(0,360,30),frac=1.05 #position of the y ticks 
        )
        ax.set_theta_direction(-1)
        
        # project strike distribution as histogram bars
           
        bars = ax.bar(theta,countbis, width=width,bottom=0.0)
        
        for r,bar in zip(countbis, bars):
            bar.set_facecolor(cm.jet(r/10.))
            bar.set_alpha(0.5)
       
        ax.set_title(Station+'_'+Component, va='bottom', fontsize=9,fontweight='bold'  # some space below the title
             )
        if st =='1':   
            ax.set_ylim(0,30)
            ax.set_yticks(np.arange(0.000,20,5)) 
            ax.set_yticklabels(np.arange(0.000,20,5),fontsize = 7)
        elif st =='3':
            ax.set_ylim(0,20)
            ax.set_yticks(np.arange(0.000,20,5)) 
            ax.set_yticklabels(np.arange(0.000,20,5),fontsize = 7)
        elif st =='7':
            ax.set_ylim(0,6)
            ax.set_yticks(np.arange(0.000,6,1)) 
            ax.set_yticklabels(np.arange(0.000,6,1),fontsize = 7)
        else :
            ax.set_ylim(0,10)
            ax.set_yticks(np.arange(0.000,10,2)) 
            ax.set_yticklabels(np.arange(0.000,10,2),fontsize = 7)
#        print 'countbis', countbis,theta
#        r_values = []
#        colors = []  
#        for r,bar in zip(countbis, bars):
#            r_values.append(r/float(max(countbis)))
#            #choose the last color input wich is proportional to the number of event considered~~~~~~~~~~~~~~~~
#            colors.append(cm.Blues(r_values[-1]))
#            bar.set_facecolor(colors[-1])
#        
#        # Add colorbar, make sure to specify tick locations to match desired ticklabels
#        colorlist = []
#        r_values.sort()
#        values = []
#        for val in r_values:
#            if val not in values:
#                values.append(val*float(max(countbis)))
#        
#            color = cm.Blues(val)
#            if color not in colorlist:
#                colorlist.append(color)
    plt.tight_layout()  
    plt.show()
    plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/Peak_to_Peak_Plots_Noratio/Peak_to_Peak_Plots_Noratio_Correlation/Azimuth_MeanPtP_%s_1141Hz_Ml_5_7_dist_200.png'%Station,dpi = 300)

  