# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 11:33:09 2016

@author: claire
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

DPeaktoPeakMean ={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DPGA = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DPGV = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DArias ={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DML={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DAz = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
Drdist = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}

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
for i in h: 
    k=f56[i].keys()
    for Station in k: 

        for Component in ['Z','R','T','E','N','H']:
            if f56[i][Station]['Exist_{}'.format(Component)][0]==True and f56[i][Station]['valid_{}'.format(Component)][0]==True and f56[i][Station]['Satured_{}'.format(Component)][0]==False : 
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
        for Component in ['Z','R','T','E','N','H']:
            if f67[i][Station]['Exist_{}'.format(Component)][0]==True and f67[i][Station]['valid_{}'.format(Component)][0]==True and f67[i][Station]['Satured_{}'.format(Component)][0]==False : 
                #print Component, f[i][Station]['PtP_{}'.format(Component)][:an],'sature', f[i][Station]['SNR_{}'.format(Component)][:]
                
                DML[Station][Component].append(f67[i].attrs['Ml'])
                DAz[Station][Component].append(f67[i][Station].attrs['Az'])
                Drdist[Station][Component].append(f67[i][Station].attrs['RDist'])
                DPGA[Station][Component].append(f67[i][Station]['PGA_{}'.format(Component)][0])
                DPGV[Station][Component].append(f67[i][Station]['PGV_{}'.format(Component)][0])
                DArias[Station][Component].append(f67[i][Station]['PGV_{}'.format(Component)][0])
                DPeaktoPeakMean[Station][Component].append(np.mean(np.array(f67[i][Station]['PtP_{}'.format(Component)][0:4])))
      
""" ############################################
                       ML
###############################################"""
        
fig,((ax1,ax2,ax3),(ax4,ax5,ax6))= plt.subplots(2, 3, sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.15, wspace=0.2)

ax1.scatter(DML['St1']['N'], DPeaktoPeakMean['St1']['N'], color = 'gold',marker ='+', label ='St1')
ax1.scatter(DML['St2']['N'], DPeaktoPeakMean['St2']['N'], color = 'forestgreen',marker ='+',label ='St2')
ax1.scatter(DML['St3']['N'], DPeaktoPeakMean['St3']['N'], color = 'red',marker ='+',label ='St3')
ax1.scatter(DML['St4']['N'], DPeaktoPeakMean['St4']['N'], color = 'mediumblue',marker ='+',label ='St4')
ax1.scatter(DML['St5']['N'], DPeaktoPeakMean['St5']['N'], color = 'coral',marker ='+',label ='St5')
ax1.scatter(DML['St6']['N'], DPeaktoPeakMean['St6']['N'], color = 'darkturquoise',marker ='+',label ='St6')
ax1.scatter(DML['St7']['N'], DPeaktoPeakMean['St7']['N'], color = 'k',marker ='+',label ='St7')
ax1.set_title('N')
ax1.set_ylabel("$mean\; Peak\;to\;Peak\; m.s^-1$")
   
ax2.scatter(DML['St1']['E'], DPeaktoPeakMean['St1']['E'], color = 'gold',marker ='+', label ='St1')
ax2.scatter(DML['St2']['E'], DPeaktoPeakMean['St2']['E'], color = 'forestgreen',marker ='+',label ='St2')
ax2.scatter(DML['St3']['E'], DPeaktoPeakMean['St3']['E'], color = 'red',marker ='+',label ='St3')
ax2.scatter(DML['St4']['E'], DPeaktoPeakMean['St4']['E'], color = 'mediumblue',marker ='+',label ='St4')
ax2.scatter(DML['St5']['E'], DPeaktoPeakMean['St5']['E'], color = 'coral',marker ='+',label ='St5')
ax2.scatter(DML['St6']['E'], DPeaktoPeakMean['St6']['E'], color = 'darkturquoise',marker ='+',label ='St6')
ax2.scatter(DML['St7']['E'], DPeaktoPeakMean['St7']['E'], color = 'k',marker ='+',label ='St7')
ax2.set_title('E')

ax3.scatter(DML['St1']['Z'], DPeaktoPeakMean['St1']['Z'], color = 'gold',marker ='+', label ='St1')
ax3.scatter(DML['St2']['Z'], DPeaktoPeakMean['St2']['Z'], color = 'forestgreen',marker ='+',label ='St2')
ax3.scatter(DML['St3']['Z'], DPeaktoPeakMean['St3']['Z'], color = 'red',marker ='+',label ='St3')
ax3.scatter(DML['St4']['Z'], DPeaktoPeakMean['St4']['Z'], color = 'mediumblue',marker ='+',label ='St4')
ax3.scatter(DML['St5']['Z'], DPeaktoPeakMean['St5']['Z'], color = 'coral',marker ='+',label ='St5')
ax3.scatter(DML['St6']['Z'], DPeaktoPeakMean['St6']['Z'], color = 'darkturquoise',marker ='+',label ='St6')
ax3.scatter(DML['St7']['Z'], DPeaktoPeakMean['St7']['Z'], color = 'k',marker ='+',label ='St7')
ax3.set_title('Z')
ax3.set_ylim(0,0.03)
ax4.scatter(DML['St1']['R'], DPeaktoPeakMean['St1']['R'], color = 'gold',marker ='+', label ='St1')
ax4.scatter(DML['St2']['R'], DPeaktoPeakMean['St2']['R'], color = 'forestgreen',marker ='+',label ='St2')
ax4.scatter(DML['St3']['R'], DPeaktoPeakMean['St3']['R'], color = 'red',marker ='+',label ='St3')
ax4.scatter(DML['St4']['R'], DPeaktoPeakMean['St4']['R'], color = 'mediumblue',marker ='+',label ='St4')
ax4.scatter(DML['St5']['R'], DPeaktoPeakMean['St5']['R'], color = 'coral',marker ='+',label ='St5')
ax4.scatter(DML['St6']['R'], DPeaktoPeakMean['St6']['R'], color = 'darkturquoise',marker ='+',label ='St6')
ax4.scatter(DML['St7']['R'], DPeaktoPeakMean['St7']['R'], color = 'k',marker ='+',label ='St7')
ax4.set_title('R')
ax4.set_ylabel("$mean\;  Peak\;to\;Peak\; m.s^-1$")
ax4.set_xlabel('$ML$')

ax5.scatter(DML['St1']['T'], DPeaktoPeakMean['St1']['T'], color = 'gold',marker ='+', label ='St1')
ax5.scatter(DML['St2']['T'], DPeaktoPeakMean['St2']['T'], color = 'forestgreen',marker ='+',label ='St2')
ax5.scatter(DML['St3']['T'], DPeaktoPeakMean['St3']['T'], color = 'red',marker ='+',label ='St3')
ax5.scatter(DML['St4']['T'], DPeaktoPeakMean['St4']['T'], color = 'mediumblue',marker ='+',label ='St4')
ax5.scatter(DML['St5']['T'], DPeaktoPeakMean['St5']['T'], color = 'coral',marker ='+',label ='St5')
ax5.scatter(DML['St6']['T'], DPeaktoPeakMean['St6']['T'], color = 'darkturquoise',marker ='+',label ='St6')
ax5.scatter(DML['St7']['T'], DPeaktoPeakMean['St7']['T'], color = 'k',marker ='+',label ='St7')
ax5.set_title('T')
ax5.set_xlabel('$ML$')

ax6.scatter(DML['St1']['H'], DPeaktoPeakMean['St1']['H'], color = 'gold',marker ='+', label ='St1')
ax6.scatter(DML['St2']['H'], DPeaktoPeakMean['St2']['H'], color = 'forestgreen',marker ='+',label ='St2')
ax6.scatter(DML['St3']['H'], DPeaktoPeakMean['St3']['H'], color = 'red',marker ='+',label ='St3')
ax6.scatter(DML['St4']['H'], DPeaktoPeakMean['St4']['H'], color = 'mediumblue',marker ='+',label ='St4')
ax6.scatter(DML['St5']['H'], DPeaktoPeakMean['St5']['H'], color = 'coral',marker ='+',label ='St5')
ax6.scatter(DML['St6']['H'], DPeaktoPeakMean['St6']['H'], color = 'darkturquoise',marker ='+',label ='St6')
ax6.scatter(DML['St7']['H'], DPeaktoPeakMean['St7']['H'], color = 'k',marker ='+',label ='St7')
ax6.set_title('H')
ax6.set_xlabel('ML')

ax6.legend(loc=(1.05, 0.7),fontsize=10)
fig.subplots_adjust(right=1.8)
plt.ylim(0,0.03)
plt.tight_layout()
plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/Peak_to_Peak_Plots_Noratio/Peak_to_Peak_Plots_Noratio_Correlation/ML_Mean_PtP_Ml_5_7Dist_200.png',bbox_inches='tight', dpi =300)

""" ############################################
                        PGA
###############################################"""
                        
fig,((ax1,ax2,ax3),(ax4,ax5,ax6))= plt.subplots(2, 3, sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.15, wspace=0.1)

ax1.scatter(DPGA['St1']['N'], DPeaktoPeakMean['St1']['N'], color = 'gold',marker ='+', label ='St1')
ax1.scatter(DPGA['St2']['N'], DPeaktoPeakMean['St2']['N'], color = 'forestgreen',marker ='+',label ='St2')
ax1.scatter(DPGA['St3']['N'], DPeaktoPeakMean['St3']['N'], color = 'red',marker ='+',label ='St3')
ax1.scatter(DPGA['St4']['N'], DPeaktoPeakMean['St4']['N'], color = 'mediumblue',marker ='+',label ='St4')
ax1.scatter(DPGA['St5']['N'], DPeaktoPeakMean['St5']['N'], color = 'coral',marker ='+',label ='St5')
ax1.scatter(DPGA['St6']['N'], DPeaktoPeakMean['St6']['N'], color = 'darkturquoise',marker ='+',label ='St6')
ax1.scatter(DPGA['St7']['N'], DPeaktoPeakMean['St7']['N'], color = 'k',marker ='+',label ='St7')
ax1.set_title('N')
ax1.set_ylabel("$mean\; Peak\;to\;Peak\; m.s^-1$")
ax1.set_xlim(0,0.015)   
ax2.scatter(DPGA['St1']['E'], DPeaktoPeakMean['St1']['E'], color = 'gold',marker ='+', label ='St1')
ax2.scatter(DPGA['St2']['E'], DPeaktoPeakMean['St2']['E'], color = 'forestgreen',marker ='+',label ='St2')
ax2.scatter(DPGA['St3']['E'], DPeaktoPeakMean['St3']['E'], color = 'red',marker ='+',label ='St3')
ax2.scatter(DPGA['St4']['E'], DPeaktoPeakMean['St4']['E'], color = 'mediumblue',marker ='+',label ='St4')
ax2.scatter(DPGA['St5']['E'], DPeaktoPeakMean['St5']['E'], color = 'coral',marker ='+',label ='St5')
ax2.scatter(DPGA['St6']['E'], DPeaktoPeakMean['St6']['E'], color = 'darkturquoise',marker ='+',label ='St6')
ax2.scatter(DPGA['St7']['E'], DPeaktoPeakMean['St7']['E'], color = 'k',marker ='+',label ='St7')
ax2.set_title('E')
ax2.set_xlim(0,0.015)
ax3.scatter(DPGA['St1']['Z'], DPeaktoPeakMean['St1']['Z'], color = 'gold',marker ='+', label ='St1')
ax3.scatter(DPGA['St2']['Z'], DPeaktoPeakMean['St2']['Z'], color = 'forestgreen',marker ='+',label ='St2')
ax3.scatter(DPGA['St3']['Z'], DPeaktoPeakMean['St3']['Z'], color = 'red',marker ='+',label ='St3')
ax3.scatter(DPGA['St4']['Z'], DPeaktoPeakMean['St4']['Z'], color = 'mediumblue',marker ='+',label ='St4')
ax3.scatter(DPGA['St5']['Z'], DPeaktoPeakMean['St5']['Z'], color = 'coral',marker ='+',label ='St5')
ax3.scatter(DPGA['St6']['Z'], DPeaktoPeakMean['St6']['Z'], color = 'darkturquoise',marker ='+',label ='St6')
ax3.scatter(DPGA['St7']['Z'], DPeaktoPeakMean['St7']['Z'], color = 'k',marker ='+',label ='St7')
ax3.set_title('Z')
ax3.set_ylim(0,0.03)
ax3.set_xlim(0,0.105)
ax4.scatter(DPGA['St1']['R'], DPeaktoPeakMean['St1']['R'], color = 'gold',marker ='+', label ='St1')
ax4.scatter(DPGA['St2']['R'], DPeaktoPeakMean['St2']['R'], color = 'forestgreen',marker ='+',label ='St2')
ax4.scatter(DPGA['St3']['R'], DPeaktoPeakMean['St3']['R'], color = 'red',marker ='+',label ='St3')
ax4.scatter(DPGA['St4']['R'], DPeaktoPeakMean['St4']['R'], color = 'mediumblue',marker ='+',label ='St4')
ax4.scatter(DPGA['St5']['R'], DPeaktoPeakMean['St5']['R'], color = 'coral',marker ='+',label ='St5')
ax4.scatter(DPGA['St6']['R'], DPeaktoPeakMean['St6']['R'], color = 'darkturquoise',marker ='+',label ='St6')
ax4.scatter(DPGA['St7']['R'], DPeaktoPeakMean['St7']['R'], color = 'k',marker ='+',label ='St7')
ax4.set_title('R')
ax4.set_ylabel("$mean\;  Peak\;to\;Peak\; m.s^+1$")
ax4.set_xlabel('$ PGA \; m.s^-2$')

ax5.scatter(DPGA['St1']['T'], DPeaktoPeakMean['St1']['T'], color = 'gold',marker ='+', label ='St1')
ax5.scatter(DPGA['St2']['T'], DPeaktoPeakMean['St2']['T'], color = 'forestgreen',marker ='+',label ='St2')
ax5.scatter(DPGA['St3']['T'], DPeaktoPeakMean['St3']['T'], color = 'red',marker ='+',label ='St3')
ax5.scatter(DPGA['St4']['T'], DPeaktoPeakMean['St4']['T'], color = 'mediumblue',marker ='+',label ='St4')
ax5.scatter(DPGA['St5']['T'], DPeaktoPeakMean['St5']['T'], color = 'coral',marker ='+',label ='St5')
ax5.scatter(DPGA['St6']['T'], DPeaktoPeakMean['St6']['T'], color = 'darkturquoise',marker ='+',label ='St6')
ax5.scatter(DPGA['St7']['T'], DPeaktoPeakMean['St7']['T'], color = 'k',marker ='+',label ='St7')
ax5.set_title('T')
ax5.set_xlabel('$PGA \; m.s^-12$')


ax6.scatter(DPGA['St1']['H'], DPeaktoPeakMean['St1']['H'], color = 'gold',marker ='+', label ='St1')
ax6.scatter(DPGA['St2']['H'], DPeaktoPeakMean['St2']['H'], color = 'forestgreen',marker ='+',label ='St2')
ax6.scatter(DPGA['St3']['H'], DPeaktoPeakMean['St3']['H'], color = 'red',marker ='+',label ='St3')
ax6.scatter(DPGA['St4']['H'], DPeaktoPeakMean['St4']['H'], color = 'mediumblue',marker ='+',label ='St4')
ax6.scatter(DPGA['St5']['H'], DPeaktoPeakMean['St5']['H'], color = 'coral',marker ='+',label ='St5')
ax6.scatter(DPGA['St6']['H'], DPeaktoPeakMean['St6']['H'], color = 'darkturquoise',marker ='+',label ='St6')
ax6.scatter(DPGA['St7']['H'], DPeaktoPeakMean['St7']['H'], color = 'k',marker ='+',label ='St7')
ax6.set_title('H')
ax6.set_xlabel('$PGA \; m.s^-2$')
ax6.legend(loc=(1.05, 0.7),fontsize=10)
plt.ylim(0,0.03)
plt.xlim(0,0.015)
fig.subplots_adjust(right=0.8)

plt.tight_layout()
plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/Peak_to_Peak_Plots_Noratio/Peak_to_Peak_Plots_Noratio_Correlation/PGA_Mean_PtP_Ml_5_7Dist_200.png',bbox_inches='tight', dpi =300)


"""######################################
                PGV
#########################################"""

fig,((ax1,ax2,ax3),(ax4,ax5,ax6))= plt.subplots(2, 3, sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.15, wspace=0.2)

ax1.scatter(DPGV['St1']['N'], DPeaktoPeakMean['St1']['N'], color = 'gold',marker ='+', label ='St1')
ax1.scatter(DPGV['St2']['N'], DPeaktoPeakMean['St2']['N'], color = 'forestgreen',marker ='+',label ='St2')
ax1.scatter(DPGV['St3']['N'], DPeaktoPeakMean['St3']['N'], color = 'red',marker ='+',label ='St3')
ax1.scatter(DPGV['St4']['N'], DPeaktoPeakMean['St4']['N'], color = 'mediumblue',marker ='+',label ='St4')
ax1.scatter(DPGV['St5']['N'], DPeaktoPeakMean['St5']['N'], color = 'coral',marker ='+',label ='St5')
ax1.scatter(DPGV['St6']['N'], DPeaktoPeakMean['St6']['N'], color = 'darkturquoise',marker ='+',label ='St6')
ax1.scatter(DPGV['St7']['N'], DPeaktoPeakMean['St7']['N'], color = 'k',marker ='+',label ='St7')
ax1.set_title('N')
ax1.set_ylabel("$mean\; Peak\;to\;Peak\; m.s^-1$")
ax1.set_xlim(0) 
ax1.set_ylim(0,0.03) 
ax2.scatter(DPGV['St1']['E'], DPeaktoPeakMean['St1']['E'], color = 'gold',marker ='+', label ='St1')
ax2.scatter(DPGV['St2']['E'], DPeaktoPeakMean['St2']['E'], color = 'forestgreen',marker ='+',label ='St2')
ax2.scatter(DPGV['St3']['E'], DPeaktoPeakMean['St3']['E'], color = 'red',marker ='+',label ='St3')
ax2.scatter(DPGV['St4']['E'], DPeaktoPeakMean['St4']['E'], color = 'mediumblue',marker ='+',label ='St4')
ax2.scatter(DPGV['St5']['E'], DPeaktoPeakMean['St5']['E'], color = 'coral',marker ='+',label ='St5')
ax2.scatter(DPGV['St6']['E'], DPeaktoPeakMean['St6']['E'], color = 'darkturquoise',marker ='+',label ='St6')
ax2.scatter(DPGV['St7']['E'], DPeaktoPeakMean['St7']['E'], color = 'k',marker ='+',label ='St7')
ax2.set_title('E')
ax2.set_xlim(0)
ax3.scatter(DPGV['St1']['Z'], DPeaktoPeakMean['St1']['Z'], color = 'gold',marker ='+', label ='St1')
ax3.scatter(DPGV['St2']['Z'], DPeaktoPeakMean['St2']['Z'], color = 'forestgreen',marker ='+',label ='St2')
ax3.scatter(DPGV['St3']['Z'], DPeaktoPeakMean['St3']['Z'], color = 'red',marker ='+',label ='St3')
ax3.scatter(DPGV['St4']['Z'], DPeaktoPeakMean['St4']['Z'], color = 'mediumblue',marker ='+',label ='St4')
ax3.scatter(DPGV['St5']['Z'], DPeaktoPeakMean['St5']['Z'], color = 'coral',marker ='+',label ='St5')
ax3.scatter(DPGV['St6']['Z'], DPeaktoPeakMean['St6']['Z'], color = 'darkturquoise',marker ='+',label ='St6')
ax3.scatter(DPGV['St7']['Z'], DPeaktoPeakMean['St7']['Z'], color = 'k',marker ='+',label ='St7')
ax3.set_title('Z')
ax3.set_xlim(0)
ax4.scatter(DPGV['St1']['R'], DPeaktoPeakMean['St1']['R'], color = 'gold',marker ='+', label ='St1')
ax4.scatter(DPGV['St2']['R'], DPeaktoPeakMean['St2']['R'], color = 'forestgreen',marker ='+',label ='St2')
ax4.scatter(DPGV['St3']['R'], DPeaktoPeakMean['St3']['R'], color = 'red',marker ='+',label ='St3')
ax4.scatter(DPGV['St4']['R'], DPeaktoPeakMean['St4']['R'], color = 'mediumblue',marker ='+',label ='St4')
ax4.scatter(DPGV['St5']['R'], DPeaktoPeakMean['St5']['R'], color = 'coral',marker ='+',label ='St5')
ax4.scatter(DPGV['St6']['R'], DPeaktoPeakMean['St6']['R'], color = 'darkturquoise',marker ='+',label ='St6')
ax4.scatter(DPGV['St7']['R'], DPeaktoPeakMean['St7']['R'], color = 'k',marker ='+',label ='St7')
ax4.set_title('R')
ax4.set_ylabel("$mean\;  Peak\;to\;Peak\; m.s^-1$")
ax4.set_xlabel('$ PGV \; m.s^-1$')

ax5.scatter(DPGV['St1']['T'], DPeaktoPeakMean['St1']['T'], color = 'gold',marker ='+', label ='St1')
ax5.scatter(DPGV['St2']['T'], DPeaktoPeakMean['St2']['T'], color = 'forestgreen',marker ='+',label ='St2')
ax5.scatter(DPGV['St3']['T'], DPeaktoPeakMean['St3']['T'], color = 'red',marker ='+',label ='St3')
ax5.scatter(DPGV['St4']['T'], DPeaktoPeakMean['St4']['T'], color = 'mediumblue',marker ='+',label ='St4')
ax5.scatter(DPGV['St5']['T'], DPeaktoPeakMean['St5']['T'], color = 'coral',marker ='+',label ='St5')
ax5.scatter(DPGV['St6']['T'], DPeaktoPeakMean['St6']['T'], color = 'darkturquoise',marker ='+',label ='St6')
ax5.scatter(DPGV['St7']['T'], DPeaktoPeakMean['St7']['T'], color = 'k',marker ='+',label ='St7')
ax5.set_title('T')
ax5.set_xlabel('$PGV \; m.s^-1$')


ax6.scatter(DPGV['St1']['H'], DPeaktoPeakMean['St1']['H'], color = 'gold',marker ='+', label ='St1')
ax6.scatter(DPGV['St2']['H'], DPeaktoPeakMean['St2']['H'], color = 'forestgreen',marker ='+',label ='St2')
ax6.scatter(DPGV['St3']['H'], DPeaktoPeakMean['St3']['H'], color = 'red',marker ='+',label ='St3')
ax6.scatter(DPGV['St4']['H'], DPeaktoPeakMean['St4']['H'], color = 'mediumblue',marker ='+',label ='St4')
ax6.scatter(DPGV['St5']['H'], DPeaktoPeakMean['St5']['H'], color = 'coral',marker ='+',label ='St5')
ax6.scatter(DPGV['St6']['H'], DPeaktoPeakMean['St6']['H'], color = 'darkturquoise',marker ='+',label ='St6')
ax6.scatter(DPGA['St7']['H'], DPeaktoPeakMean['St7']['H'], color = 'k',marker ='+',label ='St7')
ax6.set_title('H')
ax6.set_xlabel('$PGV \; m.s^-1$')
ax6.legend(loc=(1.05, 0.7),fontsize=10)
fig.subplots_adjust(right=0.8)
plt.tight_layout()
ax6.set_ylim(0,0.03)
ax6.set_xlim(0) 
plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/Peak_to_Peak_Plots_Noratio/Peak_to_Peak_Plots_Noratio_Correlation/PGV_Mean_PtP_Ml_5_7Dist_200.png',bbox_inches='tight', dpi =300)
""" ############################################
                        Arias
###############################################"""
fig,((ax1,ax2,ax3),(ax4,ax5,ax6))= plt.subplots(2, 3, sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.15, wspace=0.2)

ax1.scatter(DArias['St1']['N'], DPeaktoPeakMean['St1']['N'], color = 'gold',marker ='+', label ='St1')
ax1.scatter(DArias['St2']['N'], DPeaktoPeakMean['St2']['N'], color = 'forestgreen',marker ='+',label ='St2')
ax1.scatter(DArias['St3']['N'], DPeaktoPeakMean['St3']['N'], color = 'red',marker ='+',label ='St3')
ax1.scatter(DArias['St4']['N'], DPeaktoPeakMean['St4']['N'], color = 'mediumblue',marker ='+',label ='St4')
ax1.scatter(DArias['St5']['N'], DPeaktoPeakMean['St5']['N'], color = 'coral',marker ='+',label ='St5')
ax1.scatter(DArias['St6']['N'], DPeaktoPeakMean['St6']['N'], color = 'darkturquoise',marker ='+',label ='St6')
ax1.scatter(DArias['St7']['N'], DPeaktoPeakMean['St7']['N'], color = 'k',marker ='+',label ='St7')
ax1.set_title('N')
ax1.set_ylabel("$mean\; Peak\;to\;Peak\; m.s^-1$")
ax1.set_xlim(0,0.03) 
ax1.set_ylim(0)     
ax2.scatter(DArias['St1']['E'], DPeaktoPeakMean['St1']['E'], color = 'gold',marker ='+', label ='St1')
ax2.scatter(DArias['St2']['E'], DPeaktoPeakMean['St2']['E'], color = 'forestgreen',marker ='+',label ='St2')
ax2.scatter(DArias['St3']['E'], DPeaktoPeakMean['St3']['E'], color = 'red',marker ='+',label ='St3')
ax2.scatter(DArias['St4']['E'], DPeaktoPeakMean['St4']['E'], color = 'mediumblue',marker ='+',label ='St4')
ax2.scatter(DArias['St5']['E'], DPeaktoPeakMean['St5']['E'], color = 'coral',marker ='+',label ='St5')
ax2.scatter(DArias['St6']['E'], DPeaktoPeakMean['St6']['E'], color = 'darkturquoise',marker ='+',label ='St6')
ax2.scatter(DArias['St7']['E'], DPeaktoPeakMean['St7']['E'], color = 'k',marker ='+',label ='St7')
ax2.set_title('E')
ax2.set_xlim(0) 

ax3.scatter(DArias['St1']['Z'], DPeaktoPeakMean['St1']['Z'], color = 'gold',marker ='+', label ='St1')
ax3.scatter(DArias['St2']['Z'], DPeaktoPeakMean['St2']['Z'], color = 'forestgreen',marker ='+',label ='St2')
ax3.scatter(DArias['St3']['Z'], DPeaktoPeakMean['St3']['Z'], color = 'red',marker ='+',label ='St3')
ax3.scatter(DArias['St4']['Z'], DPeaktoPeakMean['St4']['Z'], color = 'mediumblue',marker ='+',label ='St4')
ax3.scatter(DArias['St5']['Z'], DPeaktoPeakMean['St5']['Z'], color = 'coral',marker ='+',label ='St5')
ax3.scatter(DArias['St6']['Z'], DPeaktoPeakMean['St6']['Z'], color = 'darkturquoise',marker ='+',label ='St6')
ax3.scatter(DArias['St7']['Z'], DPeaktoPeakMean['St7']['Z'], color = 'k',marker ='+',label ='St7')
ax3.set_title('Z')
ax3.set_xlim(0) 

ax4.scatter(DArias['St1']['R'], DPeaktoPeakMean['St1']['R'], color = 'gold',marker ='+', label ='St1')
ax4.scatter(DArias['St2']['R'], DPeaktoPeakMean['St2']['R'], color = 'forestgreen',marker ='+',label ='St2')
ax4.scatter(DArias['St3']['R'], DPeaktoPeakMean['St3']['R'], color = 'red',marker ='+',label ='St3')
ax4.scatter(DArias['St4']['R'], DPeaktoPeakMean['St4']['R'], color = 'mediumblue',marker ='+',label ='St4')
ax4.scatter(DArias['St5']['R'], DPeaktoPeakMean['St5']['R'], color = 'coral',marker ='+',label ='St5')
ax4.scatter(DArias['St6']['R'], DPeaktoPeakMean['St6']['R'], color = 'darkturquoise',marker ='+',label ='St6')
ax4.scatter(DArias['St7']['R'], DPeaktoPeakMean['St7']['R'], color = 'k',marker ='+',label ='St7')
ax4.set_title('R')
ax4.set_ylabel("$mean\;  Peak\;to\;Peak\; m.s^-1$")
ax4.set_xlabel('$ Arias\; m.s^-1$')

ax5.scatter(DArias['St1']['T'], DPeaktoPeakMean['St1']['T'], color = 'gold',marker ='+', label ='St1')
ax5.scatter(DArias['St2']['T'], DPeaktoPeakMean['St2']['T'], color = 'forestgreen',marker ='+',label ='St2')
ax5.scatter(DArias['St3']['T'], DPeaktoPeakMean['St3']['T'], color = 'red',marker ='+',label ='St3')
ax5.scatter(DArias['St4']['T'], DPeaktoPeakMean['St4']['T'], color = 'mediumblue',marker ='+',label ='St4')
ax5.scatter(DArias['St5']['T'], DPeaktoPeakMean['St5']['T'], color = 'coral',marker ='+',label ='St5')
ax5.scatter(DArias['St6']['T'], DPeaktoPeakMean['St6']['T'], color = 'darkturquoise',marker ='+',label ='St6')
ax5.scatter(DArias['St7']['T'], DPeaktoPeakMean['St7']['T'], color = 'k',marker ='+',label ='St7')
ax5.set_title('T')
ax5.set_xlabel('$Arias\; m.s^-1$')


ax6.scatter(DArias['St1']['H'], DPeaktoPeakMean['St1']['H'], color = 'gold',marker ='+', label ='St1')
ax6.scatter(DArias['St2']['H'], DPeaktoPeakMean['St2']['H'], color = 'forestgreen',marker ='+',label ='St2')
ax6.scatter(DArias['St3']['H'], DPeaktoPeakMean['St3']['H'], color = 'red',marker ='+',label ='St3')
ax6.scatter(DArias['St4']['H'], DPeaktoPeakMean['St4']['H'], color = 'mediumblue',marker ='+',label ='St4')
ax6.scatter(DArias['St5']['H'], DPeaktoPeakMean['St5']['H'], color = 'coral',marker ='+',label ='St5')
ax6.scatter(DArias['St6']['H'], DPeaktoPeakMean['St6']['H'], color = 'darkturquoise',marker ='+',label ='St6')
ax6.scatter(DArias['St7']['H'], DPeaktoPeakMean['St7']['H'], color = 'k',marker ='+',label ='St7')
ax6.set_title('H')
ax6.set_xlabel('$Arias\; m.s^-1$')
ax6.legend(loc=(1.05, 0.7),fontsize=10)
fig.subplots_adjust(right=0.8)
plt.tight_layout()
ax6.set_ylim(0,0.03)
ax6.set_xlim(0)
plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/Peak_to_Peak_Plots_Noratio/Peak_to_Peak_Plots_Noratio_Correlation/Arias_Mean_PtP_Ml_5_7Dist_200.png',bbox_inches='tight', dpi =300)        

""" ############################################
                        rdist
###############################################"""
                        
fig,((ax1,ax2,ax3),(ax4,ax5,ax6))= plt.subplots(2, 3, sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.15, wspace=0.2)

ax1.scatter(Drdist['St1']['N'], DPeaktoPeakMean['St1']['N'], color = 'gold',marker ='+', label ='St1')
ax1.scatter(Drdist['St2']['N'], DPeaktoPeakMean['St2']['N'], color = 'forestgreen',marker ='+',label ='St2')
ax1.scatter(Drdist['St3']['N'], DPeaktoPeakMean['St3']['N'], color = 'red',marker ='+',label ='St3')
ax1.scatter(Drdist['St4']['N'], DPeaktoPeakMean['St4']['N'], color = 'mediumblue',marker ='+',label ='St4')
ax1.scatter(Drdist['St5']['N'], DPeaktoPeakMean['St5']['N'], color = 'coral',marker ='+',label ='St5')
ax1.scatter(Drdist['St6']['N'], DPeaktoPeakMean['St6']['N'], color = 'darkturquoise',marker ='+',label ='St6')
ax1.scatter(Drdist['St7']['N'], DPeaktoPeakMean['St7']['N'], color = 'k',marker ='+',label ='St7')
ax1.set_title('N')
ax1.set_ylabel("$mean\; Peak\;to\;Peak\; m.s^-1$")
ax1.set_xlim(0) 
ax1.set_ylim(0,0.03)     
ax2.scatter(Drdist['St1']['E'], DPeaktoPeakMean['St1']['E'], color = 'gold',marker ='+', label ='St1')
ax2.scatter(Drdist['St2']['E'], DPeaktoPeakMean['St2']['E'], color = 'forestgreen',marker ='+',label ='St2')
ax2.scatter(Drdist['St3']['E'], DPeaktoPeakMean['St3']['E'], color = 'red',marker ='+',label ='St3')
ax2.scatter(Drdist['St4']['E'], DPeaktoPeakMean['St4']['E'], color = 'mediumblue',marker ='+',label ='St4')
ax2.scatter(Drdist['St5']['E'], DPeaktoPeakMean['St5']['E'], color = 'coral',marker ='+',label ='St5')
ax2.scatter(Drdist['St6']['E'], DPeaktoPeakMean['St6']['E'], color = 'darkturquoise',marker ='+',label ='St6')
ax2.scatter(Drdist['St7']['E'], DPeaktoPeakMean['St7']['E'], color = 'k',marker ='+',label ='St7')
ax2.set_title('E')
ax2.set_xlim(0) 

ax3.scatter(Drdist['St1']['Z'], DPeaktoPeakMean['St1']['Z'], color = 'gold',marker ='+', label ='St1')
ax3.scatter(Drdist['St2']['Z'], DPeaktoPeakMean['St2']['Z'], color = 'forestgreen',marker ='+',label ='St2')
ax3.scatter(Drdist['St3']['Z'], DPeaktoPeakMean['St3']['Z'], color = 'red',marker ='+',label ='St3')
ax3.scatter(Drdist['St4']['Z'], DPeaktoPeakMean['St4']['Z'], color = 'mediumblue',marker ='+',label ='St4')
ax3.scatter(Drdist['St5']['Z'], DPeaktoPeakMean['St5']['Z'], color = 'coral',marker ='+',label ='St5')
ax3.scatter(Drdist['St6']['Z'], DPeaktoPeakMean['St6']['Z'], color = 'darkturquoise',marker ='+',label ='St6')
ax3.scatter(Drdist['St7']['Z'], DPeaktoPeakMean['St7']['Z'], color = 'k',marker ='+',label ='St7')
ax3.set_title('Z')
ax3.set_xlim(0) 
  
ax4.scatter(Drdist['St1']['R'], DPeaktoPeakMean['St1']['R'], color = 'gold',marker ='+', label ='St1')
ax4.scatter(Drdist['St2']['R'], DPeaktoPeakMean['St2']['R'], color = 'forestgreen',marker ='+',label ='St2')
ax4.scatter(Drdist['St3']['R'], DPeaktoPeakMean['St3']['R'], color = 'red',marker ='+',label ='St3')
ax4.scatter(Drdist['St4']['R'], DPeaktoPeakMean['St4']['R'], color = 'mediumblue',marker ='+',label ='St4')
ax4.scatter(Drdist['St5']['R'], DPeaktoPeakMean['St5']['R'], color = 'coral',marker ='+',label ='St5')
ax4.scatter(Drdist['St6']['R'], DPeaktoPeakMean['St6']['R'], color = 'darkturquoise',marker ='+',label ='St6')
ax4.scatter(Drdist['St7']['R'], DPeaktoPeakMean['St7']['R'], color = 'k',marker ='+',label ='St7')
ax4.set_title('R')
ax4.set_ylabel("$mean\;  Peak\;to\;Peak\; m.s^-1$")
ax4.set_xlabel('$ rdist\; km$')

ax5.scatter(Drdist['St1']['T'], DPeaktoPeakMean['St1']['T'], color = 'gold',marker ='+', label ='St1')
ax5.scatter(Drdist['St2']['T'], DPeaktoPeakMean['St2']['T'], color = 'forestgreen',marker ='+',label ='St2')
ax5.scatter(Drdist['St3']['T'], DPeaktoPeakMean['St3']['T'], color = 'red',marker ='+',label ='St3')
ax5.scatter(Drdist['St4']['T'], DPeaktoPeakMean['St4']['T'], color = 'mediumblue',marker ='+',label ='St4')
ax5.scatter(Drdist['St5']['T'], DPeaktoPeakMean['St5']['T'], color = 'coral',marker ='+',label ='St5')
ax5.scatter(Drdist['St6']['T'], DPeaktoPeakMean['St6']['T'], color = 'darkturquoise',marker ='+',label ='St6')
ax5.scatter(Drdist['St7']['T'], DPeaktoPeakMean['St7']['T'], color = 'k',marker ='+',label ='St7')
ax5.set_title('T')
ax5.set_xlabel('$rdist\; km$')

ax6.scatter(Drdist['St1']['H'], DPeaktoPeakMean['St1']['H'], color = 'gold',marker ='+', label ='St1')
ax6.scatter(Drdist['St2']['H'], DPeaktoPeakMean['St2']['H'], color = 'forestgreen',marker ='+',label ='St2')
ax6.scatter(Drdist['St3']['H'], DPeaktoPeakMean['St3']['H'], color = 'red',marker ='+',label ='St3')
ax6.scatter(Drdist['St4']['H'], DPeaktoPeakMean['St4']['H'], color = 'mediumblue',marker ='+',label ='St4')
ax6.scatter(Drdist['St5']['H'], DPeaktoPeakMean['St5']['H'], color = 'coral',marker ='+',label ='St5')
ax6.scatter(Drdist['St6']['H'], DPeaktoPeakMean['St6']['H'], color = 'darkturquoise',marker ='+',label ='St6')
ax6.scatter(Drdist['St7']['H'], DPeaktoPeakMean['St7']['H'], color = 'k',marker ='+',label ='St7')
ax6.set_title('H')
ax6.set_xlabel('$rdist\;km$')
ax6.legend(loc=(1.05, 0.7),fontsize=10)
fig.subplots_adjust(right=0.8)
plt.tight_layout()
ax6.set_ylim(0,0.03)
ax6.set_xlim(0)
plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/Peak_to_Peak_Plots_Noratio/Peak_to_Peak_Plots_Noratio_Correlation/Rdistance_Mean_PtP_Ml_5_7Dist_200.png', bbox_inches='tight',dpi =300)
""" ############################################
                        Azimuth
###############################################"""
                        
fig,((ax1,ax2,ax3),(ax4,ax5,ax6))= plt.subplots(2, 3, sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.15, wspace=0.2)

ax1.scatter(DAz['St1']['N'], DPeaktoPeakMean['St1']['N'], color = 'gold',marker ='+', label ='St1')
ax1.scatter(DAz['St2']['N'], DPeaktoPeakMean['St2']['N'], color = 'forestgreen',marker ='+',label ='St2')
ax1.scatter(DAz['St3']['N'], DPeaktoPeakMean['St3']['N'], color = 'red',marker ='+',label ='St3')
ax1.scatter(DAz['St4']['N'], DPeaktoPeakMean['St4']['N'], color = 'mediumblue',marker ='+',label ='St4')
ax1.scatter(DAz['St5']['N'], DPeaktoPeakMean['St5']['N'], color = 'coral',marker ='+',label ='St5')
ax1.scatter(DAz['St6']['N'], DPeaktoPeakMean['St6']['N'], color = 'darkturquoise',marker ='+',label ='St6')
ax1.scatter(DAz['St7']['N'], DPeaktoPeakMean['St7']['N'], color = 'k',marker ='+',label ='St7')
ax1.set_title('N')
ax1.set_ylabel("$mean\; Peak\;to\;Peak\; m.s^-1$")
 
 
ax2.scatter(DAz['St1']['E'], DPeaktoPeakMean['St1']['E'], color = 'gold',marker ='+', label ='St1')
ax2.scatter(DAz['St2']['E'], DPeaktoPeakMean['St2']['E'], color = 'forestgreen',marker ='+',label ='St2')
ax2.scatter(DAz['St3']['E'], DPeaktoPeakMean['St3']['E'], color = 'red',marker ='+',label ='St3')
ax2.scatter(DAz['St4']['E'], DPeaktoPeakMean['St4']['E'], color = 'mediumblue',marker ='+',label ='St4')
ax2.scatter(DAz['St5']['E'], DPeaktoPeakMean['St5']['E'], color = 'coral',marker ='+',label ='St5')
ax2.scatter(DAz['St6']['E'], DPeaktoPeakMean['St6']['E'], color = 'darkturquoise',marker ='+',label ='St6')
ax2.scatter(DAz['St7']['E'], DPeaktoPeakMean['St7']['E'], color = 'k',marker ='+',label ='St7')
ax2.set_title('E')

ax3.scatter(DAz['St1']['Z'], DPeaktoPeakMean['St1']['Z'], color = 'gold',marker ='+', label ='St1')
ax3.scatter(DAz['St2']['Z'], DPeaktoPeakMean['St2']['Z'], color = 'forestgreen',marker ='+',label ='St2')
ax3.scatter(DAz['St3']['Z'], DPeaktoPeakMean['St3']['Z'], color = 'red',marker ='+',label ='St3')
ax3.scatter(DAz['St4']['Z'], DPeaktoPeakMean['St4']['Z'], color = 'mediumblue',marker ='+',label ='St4')
ax3.scatter(DAz['St5']['Z'], DPeaktoPeakMean['St5']['Z'], color = 'coral',marker ='+',label ='St5')
ax3.scatter(DAz['St6']['Z'], DPeaktoPeakMean['St6']['Z'], color = 'darkturquoise',marker ='+',label ='St6')
ax3.scatter(DAz['St7']['Z'], DPeaktoPeakMean['St7']['Z'], color = 'k',marker ='+',label ='St7')
ax3.set_title('Z')
ax3.set_ylim(0,0.03)

ax4.scatter(DAz['St1']['R'], DPeaktoPeakMean['St1']['R'], color = 'gold',marker ='+', label ='St1')
ax4.scatter(DAz['St2']['R'], DPeaktoPeakMean['St2']['R'], color = 'forestgreen',marker ='+',label ='St2')
ax4.scatter(DAz['St3']['R'], DPeaktoPeakMean['St3']['R'], color = 'red',marker ='+',label ='St3')
ax4.scatter(DAz['St4']['R'], DPeaktoPeakMean['St4']['R'], color = 'mediumblue',marker ='+',label ='St4')
ax4.scatter(DAz['St5']['R'], DPeaktoPeakMean['St5']['R'], color = 'coral',marker ='+',label ='St5')
ax4.scatter(DAz['St6']['R'], DPeaktoPeakMean['St6']['R'], color = 'darkturquoise',marker ='+',label ='St6')
ax4.scatter(DAz['St7']['R'], DPeaktoPeakMean['St7']['R'], color = 'k',marker ='+',label ='St7')
ax4.set_title('R')
ax4.set_ylabel("$mean\;  Peak\;to\;Peak\; m.s^-1$")
ax4.set_xlabel('$ Az\; deg$')
ax4.set_xlim(-190,190)
ax4.set_ylim(0,0.03)

ax5.scatter(DAz['St1']['T'], DPeaktoPeakMean['St1']['T'], color = 'gold',marker ='+', label ='St1')
ax5.scatter(DAz['St2']['T'], DPeaktoPeakMean['St2']['T'], color = 'forestgreen',marker ='+',label ='St2')
ax5.scatter(DAz['St3']['T'], DPeaktoPeakMean['St3']['T'], color = 'red',marker ='+',label ='St3')
ax5.scatter(DAz['St4']['T'], DPeaktoPeakMean['St4']['T'], color = 'mediumblue',marker ='+',label ='St4')
ax5.scatter(DAz['St5']['T'], DPeaktoPeakMean['St5']['T'], color = 'coral',marker ='+',label ='St5')
ax5.scatter(DAz['St6']['T'], DPeaktoPeakMean['St6']['T'], color = 'darkturquoise',marker ='+',label ='St6')
ax5.scatter(DAz['St7']['T'], DPeaktoPeakMean['St7']['T'], color = 'k',marker ='+',label ='St7')
ax5.set_title('T')
ax5.set_xlabel('$Az\; deg$')

ax5.set_xlim(-190,190)
ax6.scatter(DAz['St1']['H'], DPeaktoPeakMean['St1']['H'], color = 'gold',marker ='+', label ='St1')
ax6.scatter(DAz['St2']['H'], DPeaktoPeakMean['St2']['H'], color = 'forestgreen',marker ='+',label ='St2')
ax6.scatter(DAz['St3']['H'], DPeaktoPeakMean['St3']['H'], color = 'red',marker ='+',label ='St3')
ax6.scatter(DAz['St4']['H'], DPeaktoPeakMean['St4']['H'], color = 'mediumblue',marker ='+',label ='St4')
ax6.scatter(DAz['St5']['H'], DPeaktoPeakMean['St5']['H'], color = 'coral',marker ='+',label ='St5')
ax6.scatter(DAz['St6']['H'], DPeaktoPeakMean['St6']['H'], color = 'darkturquoise',marker ='+',label ='St6')
ax6.scatter(DAz['St7']['H'], DPeaktoPeakMean['St7']['H'], color = 'k',marker ='+',label ='St7')
ax6.set_title('H')
ax6.set_xlabel('$Az\;deg$')
ax6.legend(loc=(1.05, 0.7),fontsize=10)
fig.subplots_adjust(right=0.8)
plt.xlim(-190,190)
ax6.set_ylim(0)
plt.tight_layout()
plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/Peak_to_Peak_Plots_Noratio/Peak_to_Peak_Plots_Noratio_Correlation/Azimuth_Mean_PtP_Ml_5_7Dist_200.png', bbox_inches='tight',dpi =300)