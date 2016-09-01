
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 18:20:07 2016

@author: claire
"""

"""test Signal Wave picking """

import numpy as np
from obspy.core import read,Trace
from obspy.signal.filter import envelope
from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *
from Signal_Wave_Picking import *

Year='15' 
jJul = '110'
Hour = '11'
Station = '3' 
Second = 20*60
plot = True
#read event signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st = Read_event(Year,jJul,Hour,Second,Station, False)

#c= Covariance(st)
#select traces~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#w,v= eignevalues(st)
#print w, v
stcopy=st.copy()
stcopy.filter('bandpass',freqmin=0.2,freqmax=10,corners=4,zerophase=True)
trN, trE, trZ = SelectTracesNEZ(stcopy)
#Kurtosis test
T=5 #5 second moving windows 
trN = trN.trim(starttime=trN.stats.starttime + 60*20, endtime=trN.stats.starttime +60*30)

trZ= trZ.trim(starttime=trZ.stats.starttime + 60*20, endtime=trZ.stats.starttime +60*30)

trE=trE.trim(starttime=trE.stats.starttime + 60*20, endtime=trE.stats.starttime +60*30)

startsecond =20
endsecond=30
mpn, msn= WavePicking(trN,T,True)
mpz,msz = WavePicking(trZ,T,True)
mpe,mse = WavePicking(trE,T,True)

#Cf = CFkurtosis(trN,T) 
#f2 = F2(Cf)
#f3 = F3(f2)
#f3prim = MeanF3(f3,0.5)
#f4 = F4(f3)
#f4prim = F4(f3prim)
#
#minF4 = np.argmin(f4)
#min2F4 = np.argmin(f4[minF4+100:len(f4)-1])
#min3F4 = minF4+100+min2F4
#
#minF4p = np.argmin(f4prim)
#min2F4p = np.argmin(f4prim[minF4p+100:len(f4prim)-1])
#min3F4p = minF4p+100+min2F4p
#
#print min3F4,min3F4p
#j = st.select(component = "N")[0]
#plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#f, (ax0,ax1, ax2, ax3, ax4) = plt.subplots(5, sharex=True)
#x = range(len(Cf))
#print len(trN.data), len(Cf), len(x),len(f2), len(f3), len(f4)
#ax0.plot(x,trN.data)
#ax0.axvline(x=minF4,color= 'r')
#ax0.axvline(x=min3F4,color= 'orchid')
#ax0.axvline(x=min3F4p,color= 'aqua')
#ax1.plot(x,Cf)
#ax1.axvline(x=min3F4,color= 'orchid')
#ax1.axvline(x=min3F4p,color= 'aqua')
#ax2.plot(x, f2)
#ax2.axvline(x=min3F4,color= 'orchid')
#ax2.axvline(x=min3F4p,color= 'aqua')
#ax3.plot(x,f3)
#ax3.axvline(x=min3F4,color= 'orchid')
#ax3.axvline(x=min3F4p,color= 'aqua')
#ax4.plot(x, f4)
#ax4.plot(x,f4prim,'slategrey')
#ax4.axvline(x=min3F4,color= 'orchid')
#ax4.axvline(x=min3F4p,color= 'aqua')

#define if peaks are P or S waves ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~é
#DR, onsets = P_S_Characterisation(st,5,1.3, min3F4p, 0.4,startsecond, endsecond)
#DR, onsetp = P_S_Characterisation(st,5,1.3, minF4p, 0.4,startsecond, endsecond)