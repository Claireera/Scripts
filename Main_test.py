# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:20:35 2016

@author: claire
"""
from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *

 
Year, jjul, hour, Start, B= jJulian('20160206035700')
print jjul, Year, hour, 
##########################################################################
# 1. plot waveforms 
##########################################################################
#1.1 Selection of events ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Year = '16'
jJul = '036'
Hour = 19
plot = True
depth =23
#1.2 plot all day for a given station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#stN, stE,stZ = Plot_Merging_2hours_Seismograms(Year, jJul, '19', '20',  ['1','2','3','4','5','6'])
#stN = Plot_allstation('N',[Year, jJul, Hour],False,['1','2','3','4','5','6'])
#stE = Plot_allstation('E',[Year, jJul, Hour],False,['1','2','3','4','5','6'])
#stZ = Plot_allstation('Z',[Year, jJul, Hour],False,['1','2','3','4','5','6'])
#fig  =plt.figure()
#stncorrect = Stream_Correctionall(stN, False).trim(starttime = stN[0].stats.starttime + 57*60+18, endtime = stN[0].stats.starttime + 59*60+36)
#stncorrect.plot(fig=fig)
#fig  =plt.figure()
#stEcorrect = Stream_Correctionall(stE, False).trim(starttime =stE[0].stats.starttime +  57*60+18, endtime = stE[0].stats.starttime + 59*60+36)
#stEcorrect.plot(fig=fig)
#fig  =plt.figure()
#stZcorrect = Stream_Correctionall(stZ, False).trim(starttime = stZ[0].stats.starttime +57*60 +18, endtime =stZ[0].stats.starttime + 59*60+36)
#stZcorrect.plot(fig=fig)


stR, stT, stZ = PlotallRot('16','036', 60*57+36, 19, 22.940,120.6,depth,['1','2','3','4','5','6'])
#stallN = Plot_allstation('N',[Year, jJul, Hour],True,['1','2','3','4','5','6','7'])
#stallNcorrect = Stream_Correctionall(stallN, plot)
# 1.1. Plot all stations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot = input('do you want all stations trace plots')

## 3. Choose one station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Station = '1'
### 4. read event ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#st=Read_event(Year,str(jJul),10,1500, Station, False) 
### 4. instrumental correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#st_correc = Stream_Correction(st, Station, plot)
#
###5. Trim with LTA/STA method~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#st_trim=Trim(st_correc,1530,2, 0.5,10,40,False,plot)
### 1.4. Horizontal and vertical trace~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#trZ = st_trim.select(component='Z')[0]
#trH = RealtrH(st_trim,False)
#
### 1.5. Smooth Spectrum~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#frz, xz = Spectrum(trH, plot) 
#
#xzsmoothH,f = SmoothSpectrum(trH, plot)
#xzsmoothZ,f = SmoothSpectrum(trZ, plot)
#
#plt.plot(f,xzsmoothH/xzsmoothZ)
#plt.plot(frZ,SmoothSpectreH/SmoothSpectrez)
# 1.6. HVSR~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Mean HVSR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8. STD HVSR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9. plot and save figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~