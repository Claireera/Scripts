# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 11:21:05 2016

@author: claire

tets signal calclauted on the cluster 

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


MlMin= 6
MlMax = 7
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'

EQ = 1
Station = '7'
#Read File with EQ results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f = h5py.File('/home/claire/PHD/Working/Results/Peak_to_Peak/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s.hdf5'%(str(MlMin),str(MlMax),str(DistMax)),'r+')
g= f.keys()

h = copy.copy(g)
Year = f[h[EQ]].attrs['Year'] 
jJul = f[h[EQ]].attrs['JJul'] 
Hour=f[h[EQ]].attrs['Hour'] 
secondp= f[h[EQ]].attrs['Secondp'] 
seconds= f[h[EQ]].attrs['Seconds']
BAz =f[h[EQ]]['St%s'%Station].attrs['BAz']

Exist = f[h[EQ]]['St%s'%Station]['Exist_Z'][0]
print Year,jJul,Hour, 'St%s'%Station,Exist
#Read stream~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st=Read_event(Year,jJul,Hour,secondp, Station, False)

#satured ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~é
Satured = Saturation(st.select(component='Z')[0],True)

#6.2 Instrument correction~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stcorrect = Stream_Correction(st, Station, False)

#6.3 Rotation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strot = Rotation(stcorrect,BAz, False)
st = stcorrect.copy()
trN,trE,trZ = SelectTracesNEZ(st)
#P S wave peaking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~é
second_P, second_S = WavePicking2(trZ,5,secondp,seconds,True)
print  f[h[EQ]]['St%s'%Station]['SecondP_Z'][0],second_P
#6.7 Trim the signal between ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trTrim =  Trimtr(trZ,second_P,1.5, 0.65,7,15,False,True)

#6.8SNR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LimSNR = 2
valid,SNR= SNRstd(trZ, trTrim, LimSNR, second_P, True)

#Peak to Peak between [0.1,2.1]
Fmin = 0.01 + 0.2
Fmax = 2 + 2.2
trfilt=trZ.copy()
trFilt = trfilt.filter('bandpass',freqmin=Fmin,freqmax=Fmax,corners=4,zerophase=True)
PeaktoPeak2 , firstPeaktoPeak = PeaktoPeak(trFilt,1,False)
print PeaktoPeak, f[h[EQ]]['St%s'%Station]['PtP_Z'][:],f[h[EQ]]['St%s'%Station]['Freq_Z'][:]
print valid, f[h[EQ]]['St%s'%Station]['valid_Z'][0],SNR, f[h[EQ]]['St%s'%Station]['SNR_Z'][0],f[h[EQ]]['St%s'%Station]['Satured_Z'][0],Satured[0]