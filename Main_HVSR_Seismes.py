# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:00:27 2016

@author: claire
"""

"""Main HVSR on seismes """



from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *


jJultoplot = '206'
LStations = ['1','2','3','4','5','6','7']
MlMin= 3 
MlMax = 9
DistMin = 0
DistMax = 200
File = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
output =  '/home/claire/PHD/Working/Data/Event_stations_Caracteristics_Ml_%s_%s_Dist_%s_%s.txt'%(MlMin, MlMax, DistMin, DistMax)
for Station in LStations : 
    #select ligne where stations 
    EventCaractSelec, shape = SelectEventMLDist(File, Station, MlMin, MlMax, DistMin, DistMax, output)
    LHVSR = []
    outfile = '/home/claire/PHD/Working/results/HVSR/HVSR_plots/HVSR_%s.eps'%Station
    outfile1 = '/home/claire/PHD/Working/results/HVSR/HVSR_plots/HVSR_%sv.eps'%Station
    outfiletxt = '/home/claire/PHD/Working/results/HVSR/HVSR_Array/HVSR_%s.txt'%Station
    print len(EventCaractSelec), 'event caract'
    for j in xrange(len(EventCaractSelec)):
        
        #1.0 take event characteristics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Station, Year,jJul, Hour = ConvertDatestr(EventCaractSelec[j])
       Second =EventCaractSelec[j][4]
       print Second, 'second'
       
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
       print 'correct'
       #4.0 Trim signal between p theoirical arrival and end of signal defined as LTA/LTA ratio == 0.5~~~~~~~~~~
       st = Trim(stcorrect,Second,2, 0.5,10,40,False,plot)
       del stcorrect
       
       # 5.0 Absolute Vertical and Horizontal quadratic average~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       trH, trZ =  RealtrH(st,plot)
       print' H  V componement validate'
       #6.0 Smooth Spectrum: bandwidth = 110 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print "smooth V"
       frv,SV = SmoothSpec(trZ,110 , plot)
       print 'smooth H'
       frh, SH = SmoothSpec(trH, 110, plot)
       del trZ
       del trH
       print 'HVSr', type
       
       #5.0 HSVR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Hvsr = HVSR(SH,SV)
       LHVSR.append(Hvsr)
       print 'HVSR'      
    print 'lalal'
    AHVSR = np.array(LHVSR)
    # 7.0 mean and standard variation of HVSR 
    MeanSTDHVSR(AHVSR,frv, outfilev)
    MeanSTDHVSR(AHVSR,frh, outfile)
    np.savetxt(outfiletxt,AHVSR,delimiter="," )
    print Station, len(HVSR)
