# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 15:15:27 2016

@author: claire
"""

"""Main HVSR continus 
~~~~~~~~~~~~~~~~~~~~~

HVSR is computed continouly during a given period 
HVSR is the mean of the hourly HVSR calculated along the day for given hours ex [0,5][22,23]
"""
from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *
Year ='15'
LjJday = range(212,220)

#LHour = range(0,14)+range(22,24) #day light
#plageH = '6_21'
LHour = range(14,22) #night
plageH = '22_5'

jJultoplot=False
Hourtoplot = False

#LStations = ['1','2','3','4','5','6','7']
Station = '1'
outfile = '/home/claire/PHD/Working/Results/HVSR/HVSR_plots/HVSR_Continuous_%s_H_%s.eps'%(Station,plageH )
outfiletxt = '/home/claire/PHD/Working/Results/HVSR/HVSR_Arrays/HVSR_Continuous_HVSR_%s_H_%s.txt'%(Station, plageH)
outfilefrtxt = '/home/claire/PHD/Working/Results/HVSR/HVSR_Arrays/HVSR_Continuous_fr_%s_H_%s.txt'%(Station, plageH)
outfilejUltxt = '/home/claire/PHD/Working/Results/HVSR/HVSR_Arrays/HVSR_Continuous_jJul_%s_H_%s.txt'%(Station, plageH)
LStations = ['1','2','3']
for Station in LStations:

    AllHVSR = []
    Allfrh = []
    AllJJul = []
    for jJul in LjJday:
        LHVSR = []
        Lfr = []
        outfile = '/home/claire/PHD/Working/Results/HVSR/HVSR_plots/HVSR_Continuous_%s_H_%s.eps'%(Station,plageH )
        outfiletxt = '/home/claire/PHD/Working/Results/HVSR/HVSR_Arrays/HVSR_Continuous_HVSR_%s_H_%s.txt'%(Station, plageH)
        outfilefrtxt = '/home/claire/PHD/Working/Results/HVSR/HVSR_Arrays/HVSR_Continuous_fr_%s_H_%s.txt'%(Station, plageH)
        outfilejUltxt = '/home/claire/PHD/Working/Results/HVSR/HVSR_Arrays/HVSR_Continuous_jJul_%s_H_%s.txt'%(Station, plageH)
        for Hour in LHour:
            if Hour<10: 
                Hour='0'+str(Hour)
            else : 
                Hour=str(Hour)
            # 1.0 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if not  os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,str(jJul),Station,Year,str(jJul),Hour)) : 
               print 'file not existing' , jJul, Hour 
               Hvsr = np.empty((180299,0,))
               Hvsr[:]=np.NAN
               frh = Hvsr
               continue
            if jJul == jJultoplot and Hour==Hourtoplot: 
                plot =True
            else : 
                plot= False
            
            st=Read_event(Year,str(jJul),Hour,3599, Station, plot)
           
            # 2.0 instrumental correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            stcorrect = Stream_Correction(st, Station, plot)
            del st
            st = stcorrect.copy()
            del stcorrect
            
            # 3.0 Absolute Vertical and Horizontal Quadratic average component~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            trH,trZ = RealtrH(st,plot)
           
            # 4.0 Smooth Spectrum: flat :300 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print 'SMOOTH', jJul, '_',  Hour        
            frv, SV = SmoothSpectrum(trZ, plot)
            frh, SH = SmoothSpectrum(trH, plot)
         
            del trZ
            del trH
               
            # 5.0 HSVR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print 'HVSR'            
            Hvsr = HVSR(SH,SV)
            LHVSR.append(Hvsr)
            Lfr.append(frh)
        print len(Hvsr), len(frh)
        AHVSR= np.array(LHVSR)
        Afr = np.array(Lfr)
        Mean = np.nanmean(AHVSR,axis=0).tolist()
        fr = np.nanmean(Afr, axis=0).tolist()
        
        AllHVSR.append(Mean)
        
        Allfrh.append(fr)
        AllJJul.append(jJul)
   
        # 6.0 Save HVSR means fr and Julain day arrays for a given station 
    AllHVSR = np.array([i for i in AllHVSR]).T    
    np.savetxt(outfiletxt,AllHVSR,delimiter="," )
    np.savetxt(outfilefrtxt,Allfrh, delimiter =",")
    np.savetxt(outfilejUltxt, AllJJul, delimiter =",")
    PlotSpectrogram(Allfrh[0], AllHVSR, AllJJul, outfile,colormap="jet")
Allfrh = np.loadtxt(outfilefrtxt, delimiter =',')
AllHVSR = np.loadtxt(outfiletxt, delimiter =',')
AllJJul=  np.loadtxt(outfilejUltxt,  delimiter =',')