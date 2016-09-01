# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 20:27:42 2016

@author: claire
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 16:56:39 2016

@author: claire

Plot Peak to Peak Mean value between 1.1 to 4.1Hz for the station 3 
"""
import h5py
import numpy as np

from MAP_esults import *
import copy
lonStart= [119.5]
lonEnd =[123.3]
latStart =[21.5]
latEnd = [25.5]
name = ['Taiwan']
Long56 = [119.7,120,122.5]
Lat56 = [22.5,23,23.7]
outfile = 'MAP_EQ_4_7distmax_200'

outfile3 = 'MAP_Mean_peaktoPeak_St3_EQ_5_7distmax_200'
Dcolor = {'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
DCOLORS = {'3':'gold','3.5':'darkgoldenrod','4' :'red','4.5':'darkred' ,'5': 'turquoise','5.5':'mediumblue' ,'6':'forestgreen' ,'6.5':'k'}

DCOLORS = {'3':'turquoise','3.5':'mediumblue','4' :'limegreen','4.5':'forestgreen' ,'5':'gold' ,'5.5': 'darkgoldenrod','6':'red' ,'6.5':'darkred'}
LStations = ['1','2','3','4','5','6','7']

Long=[]
Lat=[]
MLcolor=[]
Lvalues = []
Lsize =[]
"""######################
            ML 5-6
#########################"""
MlMin= 5
MlMax = 6
DistMin = 0
DistMax = 200

f56 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s.hdf5'%(str(MlMin),str(MlMax),str(DistMax)),'r+')
g= f56.keys()
h=copy.copy(g)
for i in g:
    k=f56[i].keys()

    
    for st in ['St3','St5','St1'] :
        #print st, len(k)
        for Component in ['Z','E','N'] : 
            
            if f56[i][st]['Exist_{}'.format(Component)][0]==False or f56[i][st]['valid_{}'.format(Component)][0]==False or f56[i][st]['Satured_{}'.format(Component)][0]==True:
                #print 'siganl satured or SNR to low', i , Component,st
                if i in h : 
                    h.remove(i)
for i in h:
    k=f56[i].keys()
    if len(k)<>7 :
        print i, k,"lll"
    
    else:
        if f56[i].attrs['Ml']-int(f56[i].attrs['Ml'])>0.5 : 
           col=str(int(f56[i].attrs['Ml'])+0.5)
        else :
            col=str(int(f56[i].attrs['Ml']))
        MLcolor.append(DCOLORS[col])
        
        a = np.array(f56[i]['St3']['PtP_{}'.format('E')][0:3])/np.array(f56[i]['St5']['PtP_{}'.format('E')][0:3])
        if isnan(a[0])==False :
            #size = int(np.nanmean(np.array(f56[i]['St3']['PtP_{}'.format('E')][0:4])/(np.array(f56[i]['St5']['PtP_{}'.format('E')][0:4])+np.array(f56[i]['St1']['PtP_{}'.format('E')][0:4]))/2))*20+40
            size = int(np.nanmean(np.array(f56[i]['St3']['PtP_{}'.format('E')][0:4])))*20+40
            Lsize.append(size)
            #Lvalues.append(np.nanmean(np.array(f56[i]['St3']['PtP_{}'.format('E')][0:4])/(np.array(f56[i]['St5']['PtP_{}'.format('E')][0:4]+np.array(f56[i]['St1']['PtP_{}'.format('E')][0:4]))/2)))
            Lvalues.append(np.nanmean(np.array(f56[i]['St3']['PtP_{}'.format('E')][0:4])))
            Long.append(f56[i].attrs['Long'])
      
            Lat.append(f56[i].attrs['Lat'])
        else :
            continue
        
"""######################
            ML 6-7
#########################"""
MlMin= 6
MlMax = 7
DistMin = 0
DistMax = 200

f67 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s.hdf5'%(str(MlMin),str(MlMax),str(DistMax)),'r+')
g= f67.keys()
h=copy.copy(g)
for i in g:
    k=f67[i].keys()

    for st in ['St3','St5','St1'] :
        #print st, len(k)
        for Component in ['Z','E','N'] : 
            
            if f67[i][st]['Exist_{}'.format(Component)][0]==False or f67[i][st]['Satured_{}'.format(Component)][0]==True or f67[i][st]['valid_{}'.format(Component)][0]==False :
                print 'siganl satured or SNR to low', i , Component,st
                if i in h : 
  
                  h.remove(i)
                  
print h,'hhhhhh'                  
for i in h:
    k=f67[i].keys()
    if len(k)<>7 :
        print i, k,"lll"
    
    else:
        if f67[i].attrs['Ml']-int(f67[i].attrs['Ml'])>0.5 : 
           col=str(int(f67[i].attrs['Ml'])+0.5)
        else :
            col=str(int(f67[i].attrs['Ml']))
        MLcolor.append(DCOLORS[col])
        #size = int(np.nanmean(np.array(f67[i]['St3']['PtP_{}'.format('E')][0:4])/(np.array(f67[i]['St5']['PtP_{}'.format('E')][0:4])+np.array(f67[i]['St1']['PtP_{}'.format('E')][0:4]))/2))*10+40
        size = int(np.nanmean(np.array(f67[i]['St3']['PtP_{}'.format('E')][0:4])))*10+40
        Lsize.append(size)
        #Lvalues.append(np.nanmean(np.array(f67[i]['St3']['PtP_{}'.format('E')][0:4])/(np.array(f67[i]['St5']['PtP_{}'.format('E')][0:4])+np.array(f67[i]['St1']['PtP_{}'.format('E')][0:4]))/2))
        Lvalues.append(np.nanmean(np.array(f67[i]['St3']['PtP_{}'.format('E')][0:4])))
        Long.append(f67[i].attrs['Long'])
        Lat.append(f67[i].attrs['Lat'])
                    
#"""######################
#            ML 4-5
##########################"""
#MlMin= 4
#MlMax = 5
#DistMin = 0
#DistMax = 200
#
#f45 = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V2_Ml_%s_%s_distmax_%s.hdf5'%(str(MlMin),str(MlMax),str(DistMax)),'r+')
#g= f45.keys()
#h=copy.copy(g)
#for i in g:
#   
#   
#    for st in k :
#        print st, len(k)
#        for Component in ['Z','E','N'] : 
#            
#            if f45[i][st]['SNR_{}'.format(Component)][0]<1 :
#                print 'siganl satured or SNR to low', i , Component,st
#                if i in h : 
#                    h.remove(i)
#for i in h:
#    k=f45[i].keys()
#    if len(k)<>7 :
#        print i, k,"lll"
#    
#    else:
#        if f45[i].attrs['Ml']-int(f45[i].attrs['Ml'])>0.5 : 
#            col=str(int(f45[i].attrs['Ml'])+0.5)
#        else :
#            col=str(int(f45[i].attrs['Ml']))
#        MLcolor.append(DCOLORS[col])
#        Lvalues.append(np.mean(np.array(f45[i]['St3']['PtP_{}'.format('E')][0:3])/np.array(f45[i]['St5']['PtP_{}'.format('E')][0:3])))
#        print np.array(f45[i]['St3']['PtP_{}'.format('E')][0:3])
#        size = np.mean(np.array(f45[i]['St3']['PtP_{}'.format('E')][0:3])/np.array(f45[i]['St5']['PtP_{}'.format('E')][0:3]))*10+10
#        Lsize.append(size)
#        Long.append(f45[i].attrs['Long'])
#        Lat.append(f45[i].attrs['Lat'])
#makeMap(lonStart[0],lonEnd[0],latStart[0],latEnd[0],name[0],Long56,Lat56,outfile)#
legend='PtP'
makeMapValues(lonStart[0],lonEnd[0],latStart[0],latEnd[0],name[0],Long,Lat,Lvalues,Lsize,legend,outfile3)