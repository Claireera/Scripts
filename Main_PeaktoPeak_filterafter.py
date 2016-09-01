# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 07:35:20 2016

@author: claire
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 19:55:59 2016

@author: claire
"""

"""Main peak to peak"""

from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *


jJultoplot = '365'
LStations = ['1','2','3','4','5','6','7']

Dcolor={'1':'#2ba70f','2':'#3342ff','3':'#ff3333','4':'#ffac33',,'5':'#FF0040','6':'#4B0082','7':'#3A0235'}
MlMin= 3 
MlMax = 9
DistMin = 0
DistMax = 200
Fmin = 3
Fmax = 20
File = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
output =  '/home/claire/PHD/Working/Data/Event_stations_Caracteristics_Ml_%s_%s_Dist_%s_%s.txt'%(MlMin, MlMax, DistMin, DistMax)


for Station in LStations : 
    #select ligne where stations 
    EventCaractSelec, shape = SelectEventMLDist(File, Station, MlMin, MlMax, DistMin, DistMax, output)
    LHVSR = []
    outfile = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_plots/PeaktoPeak_%s_PBbefQsum.eps'%Station
    outfiletxt = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBbefQsum_Lptp.txt'%Station
    outfiletxt1 = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBbefQsum_LjJUl.txt'%Station
    Lptp = []
    LJJul =[]
    for j in xrange(len(EventCaractSelec)):
        
        #1.0 take event characteristics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Station, Year,jJul, Hour , Second = ConvertDatestr(EventCaractSelec[j])
       BAz = EventCaractSelec[j][8]
       print Second, 'second', BAz, 'backaz'
       
       # 2.0 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print jJul,  'read'      
       if not  os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
           print 'file not existing' , jJul, Hour 
           continue
       
       if jJul == jJultoplot: 
           plot =True
       else : 
           plot= False
       st=Read_event(Year,jJul,Hour,Second, Station, plot)
       
       #3.0 instrumental correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       stcorrect = Stream_Correction(st, Station, False)
       del st
      
       #4.0 Trim signal between p theoirical arrival and end of signal defined as LTA/LTA ratio == 0.5~~~~~~~~~~
       st = Trim(stcorrect,Second,2, 0.5,10,40,False,plot)
       del stcorrect
       
       # 5.0 Rotation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print 'rotation'       
       strot = Rotation(st,BAz, plot)
       
       # 6.0 Passband filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       strot = strot.filter('bandpass',freqmin=Fmin,freqmax=Fmax,corners=4,zerophase=True)
       
       # 7.0 Quadratic sun ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print 'rotation' 
       tSumtrace = RMS(strot, plot)
       
       
       # 7.0 Peak to Peak ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print 'ptp'
       ptpmax,ptpFirst= PeaktoPeak(tSumtracefilt, 1, plot)
       print ptpFirst, 'ptp'
       Lptp.append(ptpmax)
       LJJul.append([jJul, Hour])
    #8.0 Save arrays~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    np.savetxt(outfiletxt, np.array(Lptp))
    np.savetxt(outfiletxt1,np.array(LJJul)[:,0],fmt='%s')
    #9.0 Plot
    plt.plot(np.array(LJJul)[:,0],Lptp, marker = '.', markersize = '2', color =Dcolor[Station], label= str(Station))


plt.ylabel('Peak to Peak')
plt.xlabel('Julian day')
plt.legend()
plt.savefig(outfile)

