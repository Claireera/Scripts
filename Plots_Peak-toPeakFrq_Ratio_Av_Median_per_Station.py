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
DSNR={'St1':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St2':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St3':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St4':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St5':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St6':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]},'St7':{'N':[],'E':[],'R':[],'T':[],'H':[],'Z':[]}}
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


print 'Dictionnary done' 

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
frqmin = [k/10. for k in range(2,110,10)]
frqmax = range(2,13,1)
frqcentral = [(frqmin[i]+frqmax[i])/2 for i in xrange(len(frqmax))]
## Setting up a colormap that's a simple transtion

#median~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
fig,  axes = plt.subplots(2, 3, sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.15, wspace=0.05)
#dtN = pd.DataFrame(np.array([np.median(DPeaktoPeak_Ratio['St1']['N'],axis=0),
#                            np.median(DPeaktoPeak_Ratio['St2']['N'],axis=0),
#                            np.median(DPeaktoPeak_Ratio['St3']['N'],axis=0),
#                            np.median(DPeaktoPeak_Ratio['St4']['N'],axis=0),
#                            np.median(DPeaktoPeak_Ratio['St5']['N'],axis=0),
#                            np.median(DPeaktoPeak_Ratio['St6']['N'],axis=0),
#                            np.median(DPeaktoPeak_Ratio['St7']['N'],axis=0)]).T,index = frqcentral,columns = ['St1','St2','St3','St4','St5','St6','St7'])

dtN = pd.DataFrame(np.array([np.median(DPeaktoPeak_Ratio['St1']['N'],axis=0),
                            np.median(DPeaktoPeak_Ratio['St2']['N'],axis=0),
                            np.median(DPeaktoPeak_Ratio['St3']['N'],axis=0),
                            np.median(DPeaktoPeak_Ratio['St5']['N'],axis=0),
                            np.median(DPeaktoPeak_Ratio['St6']['N'],axis=0),
                            ]).T,index = frqcentral,columns = ['St1','St2','St3','St5','St6'])

ax1 = dtN.plot(ax=axes[0,0], title = '$N$',legend=False,color=['gold','forestgreen' ,'red' , 'mediumblue','darkturquoise'])
ax1.set_ylabel("$median\; Peak\;to\;Peak\; ratio_{5}\;m.s^-1$")
   
#dtE = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['E'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St2']['E'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St3']['E'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St4']['E'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St5']['E'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St6']['E'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St7']['E'],axis=0)]).T,index = frqcentral,columns = ['St1','St2','St3','St4','St5','St6','St7'])
                            
dtE = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['E'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St2']['E'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St3']['E'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St5']['E'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St6']['E'],axis=0),
                            ]).T,index = frqcentral,columns = ['St1','St2','St3','St5','St6'])
dtE.plot(ax=axes[0,1], title = '$E$',legend=False,color=['gold','forestgreen' ,'red' , 'mediumblue','darkturquoise'])
 
#dtZ = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['Z'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St2']['Z'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St3']['Z'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St4']['Z'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St5']['Z'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St6']['Z'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St7']['Z'],axis=0)]).T,index = frqcentral,columns = ['St1','St2','St3','St4','St5','St6','St7'])
                            
dtZ = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['Z'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St2']['Z'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St3']['Z'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St5']['Z'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St6']['Z'],axis=0),
                            ]).T,index = frqcentral,columns = ['St1','St2','St3','St5','St6'])
dtZ.plot(ax=axes[0,2],  title = '$Z$',legend=False,color=['gold','forestgreen' ,'red' , 'mediumblue','darkturquoise'])
 
#dtR = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['R'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St2']['R'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St3']['R'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St4']['R'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St5']['R'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St6']['R'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St7']['R'],axis=0)]).T,index = frqcentral,columns = ['St1','St2','St3','St4','St5','St6','St7'])
                            
dtR = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['R'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St2']['R'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St3']['R'],axis=0),
                            
                            np.median( DPeaktoPeak_Ratio['St5']['R'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St6']['R'],axis=0),
                           ]).T,index = frqcentral,columns = ['St1','St2','St3','St5','St6'])
ax2 = dtR.plot(ax=axes[1,0], title = '$R$',legend=False,color=['gold','forestgreen' ,'red' , 'mediumblue','darkturquoise'])
ax2.set_ylabel("$median\; Peak\;to\;Peak\; ratio_{5}\; m.s^-1$")
ax2.set_xlabel('$Frequencies\; Hz$')

#dtT = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['T'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St2']['T'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St3']['T'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St4']['T'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St5']['T'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St6']['T'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St7']['T'],axis=0)]).T,index = frqcentral,columns = ['St1','St2','St3','St4','St5','St6','St7'])

dtT = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['T'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St2']['T'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St3']['T'],axis=0),
                         
                            np.median( DPeaktoPeak_Ratio['St5']['T'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St6']['T'],axis=0),]).T,index = frqcentral,columns = ['St1','St2','St3','St5','St6'])
ax3 = dtT.plot(ax=axes[1,1], title = '$T$',legend=False,color=['gold','forestgreen' ,'red' , 'mediumblue','darkturquoise'])
ax3.set_xlabel('$Frequencies\; Hz$')

#dtH = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['H'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St2']['H'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St3']['H'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St4']['H'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St5']['H'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St6']['H'],axis=0),
#                            np.median( DPeaktoPeak_Ratio['St7']['H'],axis=0)]).T,index = frqcentral,columns = ['$St1$','$St2$','$St3$','$St4$','$St5$','$St6$','$St7$'])

dtH = pd.DataFrame(np.array([np.median( DPeaktoPeak_Ratio['St1']['H'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St2']['H'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St3']['H'],axis=0),
                            
                            np.median( DPeaktoPeak_Ratio['St5']['H'],axis=0),
                            np.median( DPeaktoPeak_Ratio['St6']['H'],axis=0),
                            ]).T,index = frqcentral,columns = ['St1','St2','St3','St5','St6'])
 
ax4 = dtH.plot(ax=axes[1,2], title = '$H$',legend=True,color=['gold','forestgreen' ,'red' , 'mediumblue','darkturquoise'])
ax4.legend(loc=(1.05, 0.7),fontsize=10)
ax4.set_xlabel('$Frequencies\; Hz$')    
# Creating legend and title for the figure. Legend created with figlegend(), title with suptitle()
fig.subplots_adjust(right=0.8)

# If a ListedColormap is used, the length of the bounds array must be
# one greater than the length of the color list.  The bounds must be
# monotonically increasing.

plt.savefig('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/Peak_to_Peak_Plots_Ratios/PtP_RATIO_Av5_Ml4_7Dist%s_%s_median_per_Stations.png'%(DistMin,DistMax), dpi =300)     

