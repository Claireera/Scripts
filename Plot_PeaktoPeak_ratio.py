# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 09:57:47 2016

@author: claire

CAREFULL ONLY EQ WITH ALL THE REF STATIONS ARE TAKEN INTO ACCOUNT
"""

import h5py
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import matplotlib.cm  as cm 
import matplotlib as mpl
import copy

DPeaktoPeak ={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DPeaktoPeak_Ratio ={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DPGA = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DPGV = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DArias ={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DML={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
Dcolor = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DCOLORS = {'3':'turquoise','3.5':'mediumblue','4' :'limegreen','4.5':'forestgreen' ,'5':'gold' ,'5.5': 'darkgoldenrod','6':'red' ,'6.5':'darkred'}
LStations = ['1','2','3','4','5','6','7']

# VERIFICATION THAT THE SIGNAL FOR ALL THE STATIONS FOR A GIVEN eq IS VALID


"""########################################### 
            ML 5-6
##############################################""" 
MlMin= 5
MlMax = 6
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'
f56 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s.hdf5'%(str(MlMin),str(MlMax),str(DistMax)),'r+')
g= f56.keys()
h = copy.copy(g)
print 'Leng56before',len(g)
a=0
for i in g:
    a=a+1
    k=f56[i].keys()

    for st in ['St1','St3','St5','St2','St6'] :
        for Component in ['Z','E','N','R','T']  : 
            
            if f56[i][st]['Exist_{}'.format(Component)][0]==False or f56[i][st]['valid_{}'.format(Component)][0]==False or f56[i][st]['Satured_{}'.format(Component)][0]==True:
                #print 'Exist',f56[i][st]['Exist_{}'.format(Component)][0], i , Component,st,f56[i][st]['valid_{}'.format(Component)][0]
                if i in h : 
                    h.remove(i)
                        
print 'Leng56after',len(g), len(h),a
for i in h:
    for Station in ['St1','St3','St5','St2','St6'] :
        for Component in ['Z','R','T','E','N','H']:
            
            #print Component, f[i][Station]['PtP_{}'.format(Component)][:],'sature', f[i][Station]['SNR_{}'.format(Component)][:]
            DML[Station][Component].append(f56[i].attrs['Ml'])
            #print 'ML:',f56[i].attrs['Ml']
            if f56[i].attrs['Ml']-int(f56[i].attrs['Ml'])>0.5 : 
                col=str(int(f56[i].attrs['Ml'])+0.5)
            else :
                col=str(int(f56[i].attrs['Ml']))
            DML[Station][Component].append(f56[i].attrs['Ml'])
            Dcolor[Station][Component].append(DCOLORS[col]) 
            DPeaktoPeak[Station][Component].append(f56[i][Station]['PtP_{}'.format(Component)][:])
            DPGA[Station][Component].append(f56[i][Station]['PGA_{}'.format(Component)][0])
            DPGV[Station][Component].append(f56[i][Station]['PGV_{}'.format(Component)][0])
            DArias[Station][Component].append(f56[i][Station]['Arias_{}'.format(Component)][0])
        else : 
            print 'tt'
        #print i, '56 is statured or has a SNR too low', f56[i][Station]['valid_{}'.format(Component)][0],f56[i][Station]['SNR_{}'.format(Component)][0],f56[i][Station]['Satured_{}'.format(Component)][0]
f56.close() 
print len( h) , len(g), '56'

"""########################################### 
            ML 6-7
##############################################""" 
MlMin= 6
MlMax = 7
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'
f67 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s.hdf5'%(MlMin,MlMax,DistMax),'r+')
g=f67.keys()
h = copy.copy(g)
print 'Leng56before',len(g)
for i in g:
    k=f67[i].keys()

    for st in ['St1','St3','St5','St2','St6'] :
        for Component in ['Z','R','T','E','N','H'] : 
            
            if f67[i][st]['Exist_{}'.format(Component)][0]==False or f67[i][st]['valid_{}'.format(Component)][0]==False or f67[i][st]['Satured_{}'.format(Component)][0]==True:
                print 'Exist',f67[i][st]['Exist_{}'.format(Component)][0], i , Component,st,f67[i][st]['valid_{}'.format(Component)][0]
                if i in h : 
                    h.remove(i)
                        
print 'Leng67after',len(g), len(h),h 
for i in h:
    for Station in ['St1','St3','St5','St2','St6'] :
        for Component in ['Z','R','T','E','N','H']:

       
            if f67[i].attrs['Ml']-int(f67[i].attrs['Ml'])>0.5 : 
                col=str(int(f67[i].attrs['Ml'])+0.5)
            else :
                col=str(int(f67[i].attrs['Ml']))
            try :
        
                DML[Station][Component].append(f67[i].attrs['Ml'])
                Dcolor[Station][Component].append(DCOLORS[col]) 
                DPeaktoPeak[Station][Component].append(f67[i][Station]['PtP_{}'.format(Component)][:])
                DPGA[Station][Component].append(f67[i][Station]['PGA_{}'.format(Component)][0])
                DPGV[Station][Component].append(f67[i][Station]['PGV_{}'.format(Component)][0])
                DArias[Station][Component].append(f67[i][Station]['Arias_{}'.format(Component)][0])
               
            except ValueError:
                print("Oops!  %s_%s_%s was no valid File.  Try again..." %(i,Station,Component))
                continue                
        else :
            print 'tt'
            # print i, ' is statured or has a SNR too low', f67[i][Station]['valid_{}'.format(Component)][0],f67[i][Station]['SNR_{}'.format(Component)][0],f67[i][Station]['Satured_{}'.format(Component)][0]
f67.close()  

print len( h) , len(g), '67'
#"""########################################### 
#            ML 4-5
###############################################""" 
#
#MlMin= 4
#MlMax = 5
#DistMin = 0
#DistMax = 200
#jJultoplot = 206
#refStation = '5'
#f45 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V2_Ml_%s_%s_distmax_%s.hdf5'%(MlMin,MlMax,DistMax),'r+')
#g= f45.keys()
#
## VERIFICATION THAT THE SIGNAL FOR ALL THE STATIONS FOR A GIVEN eq IS VALID
#print 'len before 45', len(g)
#for i in g:
#    k=f45[i].keys()
#   
#    if len(k)<>7 :
#        print 'there are no 7 stations'
#        g.remove(i)
#    else : 
#        for st in k :
#            for Component in ['Z','R','T','E','N','H'] : 
#                if f45[i][st]['valid_{}'.format(Component)][0]==False or f45[i][st]['Satured_{}'.format(Component)][0]==True :
#                    print 'siganl satured or SNR to low', i , Component,st
#                    if i in g : 
#                        g.remove(i)
#                   
#print 'len fater45', len(g)
#
#for i in g:
#    k=f45[i].keys()
#    if len(k)==7 :
#        for Station in k: 
#            for Component in ['Z','R','T','E','N','H']:
#                                 
#                    if f45[i].attrs['Ml']-int(f45[i].attrs['Ml'])>0.5 : 
#                        col=str(int(f45[i].attrs['Ml'])+0.5)
#                    else :
#                        col=str(int(f45[i].attrs['Ml']))
#                    try: 
#                        DML[Station][Component].append(f45[i].attrs['Ml'])
#                        Dcolor[Station][Component].append(DCOLORS[col]) 
#                        DPeaktoPeak[Station][Component].append(f45[i][Station]['PtP_{}'.format(Component)][:])
#                        DPGA[Station][Component].append(f45[i][Station]['PGA_{}'.format(Component)][0])
#                        DPGV[Station][Component].append(f45[i][Station]['PGV_{}'.format(Component)][0])
#                        DArias[Station][Component].append(f45[i][Station]['Arias_{}'.format(Component)][0])
#                    except ValueError:
#                        print("Oops!  %s_%s_%s was no valid File.  Try again..." %(i,Station,Component))
#                        continue       
#             
#                # print i, ' is statured or has a SNR too low', f67[i][Station]['valid_{}'.format(Component)][0],f67[i][Station]['SNR_{}'.format(Component)][0],f67[i][Station]['Satured_{}'.format(Component)][0]
#f45.close()  
#
#"""########################################### 
#            ML 3-4
###############################################""" 
#
#MlMin= 3
#MlMax = 4
#DistMin = 0
#DistMax = 200
#jJultoplot = 206
#refStation = '5'
#f34 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V2_Ml_%s_%s_distmax_%s.hdf5'%(MlMin,MlMax,DistMax),'r+')
#g= f34.keys()
#print 'LENG BEFORE 67',len(g)
#for i in g:
#    k=f34[i].keys()
#    if len(k)<>7 :
#        print 'there are no 7 stations'
#        g.remove(i)
#    else : 
#        for st in k :
#            for Component in ['Z','R','T','E','N','H'] : 
#                if 34[i][st]['valid_{}'.format(Component)][0]==False or f34[i][st]['Satured_{}'.format(Component)][0]==True :
#                    
#                    if i in g :                     
#                        g.remove(i)
#for i in g:
#    k=f34[i].keys()
#    if len(k)==7 :
#        for Station in k: 
#            for Component in ['Z','R','T','E','N','H']:
#                                 
#                    if f34[i].attrs['Ml']-int(f34[i].attrs['Ml'])>0.5 : 
#                        col=str(int(f34[i].attrs['Ml'])+0.5)
#                    else :
#                        col=str(int(f34[i].attrs['Ml']))
#                    try: 
#                        DML[Station][Component].append(f34[i].attrs['Ml'])
#                        Dcolor[Station][Component].append(DCOLORS[col]) 
#                        DPeaktoPeak[Station][Component].append(f34[i][Station]['PtP_{}'.format(Component)][:])
#                        DPGA[Station][Component].append(f34[i][Station]['PGA_{}'.format(Component)][0])
#                        DPGV[Station][Component].append(f34[i][Station]['PGV_{}'.format(Component)][0])
#                        DArias[Station][Component].append(f34[i][Station]['Arias_{}'.format(Component)][0])
#                    except ValueError:
#                        print("Oops!  %s_%s_%s was no valid File.  Try again..." %(i,Station,Component))
#                        continue       
#           
#                # print i, ' is statured or has a SNR too low', f67[i][Station]['valid_{}'.format(Component)][0],f67[i][Station]['SNR_{}'.format(Component)][0],f67[i][Station]['Satured_{}'.format(Component)][0]
#f34.close()   

print 'Dictionnary done'
##plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frqmin = [k/10. for k in range(2,110,10)]
frqmax = range(2,13,1)
frqcentral = [(frqmin[i]+frqmax[i])/2 for i in xrange(len(frqmax))]
## Setting up a colormap that's a simple transtion
#mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
#
## Using contourf to provide my colorbar info, then clearing the figure
#Z = [[0,0],[0,0]]
#levels = range(min,max+step,step)
#CS3 = plt.contourf(Z, levels, cmap=mymap)
#plt.clf()


"""Ratio~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
DPeaktoPeak_Ratio ={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
for st in ['St2','St3','St6','St5','St1']: 
    
    for component in ['Z','R','T','E','N','H']:
        DPeaktoPeak_Ratio[st][component]= []
        for i in xrange(len(DPeaktoPeak[st][component])):
            #print i ,len(DPeaktoPeak[st][component]), component, st,len(DPeaktoPeak['St7'][component]),'7',len(DPeaktoPeak['St1'][component]),'1',len(DPeaktoPeak['St5'][component])
            Li = [DPeaktoPeak[st][component][i][j]/float(DPeaktoPeak['St1'][component][i][j] +DPeaktoPeak['St5'][component][i][j])/2 for j in xrange(len(DPeaktoPeak['St1'][component][i]))]
            #Li = [DPeaktoPeak[st][component][i][j]/float(DPeaktoPeak['St5'][component][i][j]) for j in xrange(len(DPeaktoPeak['St5'][component][i]))]
            DPeaktoPeak_Ratio[st][component].append(Li)
        
        
    
"""Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

for key in DPeaktoPeak_Ratio : 
    if key in ['St2','St3','St6','St5','St1']: 
        if len(DPeaktoPeak_Ratio[key]['N'])<>0:
            
            fig,  axes = plt.subplots(2, 3, sharex='col', sharey='row')
            fig.subplots_adjust(hspace=0.15, wspace=0.05)
            dtN = pd.DataFrame(np.array(DPeaktoPeak_Ratio[key]['N']).T,index = frqcentral)
            mycolor = Dcolor[key]['N']
            
            ax1 = dtN.plot(ax=axes[0,0], title = '${}.N$'.format(key),legend=False,color=mycolor)
            ax1.set_ylabel("$Peak\;to\;Peak\;Ratio_{15}\; m.s^-1$")
            mycolor = Dcolor[key]['E']
            dtE = pd.DataFrame(np.array(DPeaktoPeak_Ratio[key]['E']).T,index = frqcentral)
            dtE.plot(ax=axes[0,1], title = '${}.E$'.format(key),legend=False,color=mycolor)
            mycolor = Dcolor[key]['Z']
            dtZ = pd.DataFrame(np.array(DPeaktoPeak_Ratio[key]['Z']).T,index = frqcentral)
            dtZ.plot(ax=axes[0,2],  title = '${}.Z$'.format(key),legend=False,color=mycolor)
            mycolor = Dcolor[key]['R']
            dtR = pd.DataFrame(np.array(DPeaktoPeak_Ratio[key]['R']).T,index = frqcentral)
            ax2 = dtR.plot(ax=axes[1,0], title = '${}.R$'.format(key),legend=False,color=mycolor)
            ax2.set_ylabel("$Peak\;to\;Peak\;Ratio_{15}\; m.s^-1$")
            ax2.set_xlabel('$Frequencies\; Hz$')
            mycolor = Dcolor[key]['T']
            dtT = pd.DataFrame(np.array(DPeaktoPeak_Ratio[key]['T']).T,index = frqcentral)
            ax3 = dtT.plot(ax=axes[1,1], title = '${}.T$'.format(key),legend=False,color=mycolor)
            ax3.set_xlabel('$Frequencies\; Hz$')
            mycolor = Dcolor[key]['H']
            dtH = pd.DataFrame(np.array(DPeaktoPeak_Ratio[key]['H']).T,index = frqcentral)
         
            ax4 = dtH.plot(ax=axes[1,2], title = '${}.H$'.format(key),legend=False,color=mycolor)
            ax4.set_xlabel('$Frequencies\; Hz$')    
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
            cmap = mpl.colors.ListedColormap(['turquoise','mediumblue','limegreen','forestgreen' ,'gold' , 'darkgoldenrod','red' ,'darkred'])
            cmap.set_over('0.25')
            cmap.set_under('0.75')
            
            # If a ListedColormap is used, the length of the bounds array must be
            # one greater than the length of the color list.  The bounds must be
            # monotonically increasing.
            bounds = [3,3.5,4,4.5,5,5.5,6,6.5]
            norm = mpl.colors.Normalize(3,7)
            cb2 = mpl.colorbar.ColorbarBase(ax = cbar_ax,cmap=cmap,norm =norm,ticks=bounds)
      
            plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PtP_RATIO_15_AvSt5_{}_Ml5_7Dist%s_%s.png'.format(key)%(DistMin,DistMax), dpi =300)
#    #if len(k)<>7 :
#       # print 'EQ',i, 'listst',k
#        
#for j in f[g[2]]:
#    print f[g[2]][j].keys()
#    if 'PtP_N' in f[g[2]][j].keys() : 
#        print 'PtP_N' , f[g[2]][j]['PtP_N'][9]
#
#Array = np.loadtxt('/home/claire/PHD/Working/Results/A_St_EQ_ML3.1_3.1.txt')
#
#St1 = np.argwhere((Array[:,0]==1.))
#St2 = np.argwhere((Array[:,0]==2.))
#St3 = np.argwhere((Array[:,0]==3.))
#St4 = np.argwhere((Array[:,0]==4.))
#St5 = np.argwhere((Array[:,0]==5.))
#St6 = np.argwhere((Array[:,0]==6.))
#St7 = np.argwhere((Array[:,0]==7.))
