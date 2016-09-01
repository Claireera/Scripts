# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 19:54:28 2016

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
jJulref = 235

Dcolor={'1':'#2ba70f','2':'#3342ff','3':'#ff3333','4':'#ffac33','5':'#D7DF01','6':'#4B0082','7':'#FF00FF'}
MlMin= 3 
MlMax = 9
DistMin = 0
DistMax = 200
Fmin = 0.2
Fmax = 3
File = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics_TEST.txt'
output =  '/home/claire/PHD/Working/Data/Event_stations_Caracteristics_Ml_%s_%s_Dist_%s_%stets.txt'%(MlMin, MlMax, DistMin, DistMax)
#
for Station in LStations : 
    #select ligne where stations 
    EventCaractSelec, shape = SelectEventMLDistJulref(File, Station, MlMin, MlMax, DistMin, DistMax, jJulref, output)
    LHVSR = []
    outfile = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_Plots/PeaktoPeak_%s_PBafterQsum_%s%s.eps'%(Station,str(Fmin),str(Fmax))
    outfiletxtH = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUlMlPytp%s%s_H.txt'%(Station,str(Fmin),str(Fmax))
    outfiletxtR = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUlMlPytp%s%s_R.txt'%(Station,str(Fmin),str(Fmax))
    outfiletxtT = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUlMlPytp%s%s_T.txt'%(Station,str(Fmin),str(Fmax))
    outfiletxtZ = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUlMlPytp%s%s_Z.txt'%(Station,str(Fmin),str(Fmax))
    outfiletxtZabs = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUlMlPytp%s%s_Zabs.txt'%(Station,str(Fmin),str(Fmax))
    outfiletxtMaxH = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeakMax_%s_PBafterQsum_LjJUlMlPytp%s%s_H.txt'%(Station,str(Fmin),str(Fmax))
    outfiletxtMaxZabs = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeakMax_%s_PBafterQsum_LjJUlMlPytp%s%s_Zabs.txt'%(Station,str(Fmin),str(Fmax))
    Lptp = []
    LJJulPtpR =[]
    LJJulPtpT=[]
    LJJulPtpZ =[]
    LJJulPtpH =[] 
    LJJulPtpZabs =[]
    LJJulMaxH =[] 
    LJJulMaxZabs =[]
    for j in xrange(len(EventCaractSelec)):
        
        ###1.0 take event characteristics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Station, Year,jJul, Hour , Second,ml = ConvertDatestr(EventCaractSelec[j])
       BAz = EventCaractSelec[j][8]
       print  jJul, Hour , Year, 'year'
      
       ### 2.0 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print jJul,  'read'      
       if not  os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
           print 'file not existing' , jJul, Hour , Year
           continue
       
       if jJul == jJultoplot: 
           plot =True
       else : 
           plot= False
       st=Read_event(Year,jJul,Hour,Second, Station, plot)
      
       ###3.0 Instrumental correction~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       stcorrect = Stream_Correction(st, Station, plot)
       
       ###4.0 Trim signal between p theoirical arrival and end of signal defined as LTA/LTA ratio == 0.5~~~~~~~~~~
       st = Trim(stcorrect,Second,2, 0.5,10,40,False,plot)
       del stcorrect
       
       ### 5.0 Rotation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print 'rotation'       
       strot = Rotation(st,BAz, plot)
            
       ### 7.0 Passband filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       stfilt = ttrot.filter('bandpass',freqmin=Fmin,freqmax=Fmax,corners=4,zerophase=True) 
       
       ### 8.0 Trace definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
       trR, trT, trZ = SelectTracesRTZ( stfilt)
       trH, trZabs = RealtrH( stfilt,plot)
       
       ### 9.0 Peak to Peak ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print 'peak to Peak'
       #9.1 on RTZ traces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       ptpmaxZ,ptpFirstZ= PeaktoPeak(trZ, 1,plot)
       LJJulPtpZ.append([jJul, Hour, ml,ptpmaxZ])
       
       ptpmaxR,ptpFirstR= PeaktoPeak(trR, 1,plot)
       LJJulPtpR.append([jJul, Hour, ml,ptpmaxR])
       
       ptpmaxT,ptpFirstT= PeaktoPeak(trT, 1,plot)
       LJJulPtpT.append([jJul, Hour, ml,ptpmaxT])
       #9.2 on Vertical Absolute and horizontal quadratic average ~~~~~~~~~~~~~~~~~
       ptpmaxZabs,ptpFirstZabs= PeaktoPeak_positive(trZabs,plot)
       LJJulPtpZabs.append([jJul, Hour, ml,ptpmaxZabs])
       LJJulMaxZabs.append([jJul, Hour, ml,np.amax(trZabs.data)])
       
       ptpmaxH,ptpFirstH= PeaktoPeak_positive(trH, plot)
       LJJulPtpH.append([jJul, Hour, ml,ptpmaxH])
       LJJulMaxH.append([jJul, Hour, ml,np.amax(trH.data)])
       
    ###8.0 Save arrays~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # np.savetxt(outfiletxt, np.array(Lptp))
    np.savetxt(outfiletxtR,np.array(LJJulPtRp),fmt='%s')
    np.savetxt(outfiletxtZ,np.array(LJJulPtpZ),fmt='%s')
    np.savetxt(outfiletxtT,np.array(LJJulPtpT),fmt='%s')
    np.savetxt(outfiletxtZabs,np.array(LJJulPtpZabs),fmt='%s')
    
    np.savetxt(outfiletxtH,np.array(LJJulPtpH),fmt='%s')
    np.savetxt(outfiletxtMaxH,np.array(LJJulMaxH),fmt='%s')
    
    np.savetxt(outfiletxtZabs,np.array(LJJulPtpZabs),fmt='%s')
    np.savetxt(outfiletxtMaxZabs,np.array(LJJulMaxZabs),fmt='%s')
    
