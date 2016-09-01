# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 18:25:20 2016

@author: claire
"""

"""functions: 
- Read_event : return a stream containing the 3 traces N, Z,E of an event calib =1
- Stream_Correction :return a st whitening and corrected from instrumental response with calibration = 1!!!
- SelectTracesNEZ : returns horizontal traces TrE and trN  and vertical trZ
- SelectTracesRTZ : returns horizontal traces TrR and trT, and vertical TRz
- Rotation : return stream rotate thnaks to the Back azimuth between the event and the station
- Stream_PBfilter : passband filter of a stream
- RMS: return quadradic sum with all signal multiplied by their signs
- RealtrH : return quadradic sum of the horizontal components with all signal multiplied by their signs
- Stream_PBfilter : return stream filter with a passband
- LTASTA :return the cut on and off of the recursive LTA/STA list [[cuton, cutoff], [cuton, cutoff]]
- Trim : return the stream cut from the P therical arrival time, to the cut off lTA/ STA are compute on Z component for the moment 
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import read,UTCDateTime,Stream,Trace
from obspy.signal.filter import bandpass
from obspy.signal.trigger import arPick, classicSTALTA,plotTrigger,trigger_onset,recursive_sta_lta
from Events_Selections import ConvertDatestr
import sys
import math

def Read_event(Year,jJul,Hour,Second,Station, plot):
    """return a stream containing the 3 traces N, Z,E of an event calib =1
   * input :
        - station: station considered ex : '1'
        - year :  year of the event ex 15
        - jJul : julian day
        - Hour : hour of the event ex '08
        - plot : if True plot will be shown  '
    *output : 
        -st : type : stream ; Stream containing the 3 traces ENZ
    * exemples :
      1. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Station, Year,jJul, Hour = ConvertDatestr(B[j])
        ==> (Station, Year,jJul, Hour) = ('2', '15', '001', '00')
       st =  Read_event(Year,jJul,Hour,12, Station, True)
       
      2.  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       for j in [1,222344,222345]:
    
            Station, Year,jJul, Hour = ConvertDatestr(B[j])
            print ConvertDatestr(B[j])
        
           
            if not  os.path.exists('/home/claire/PHD/Working/Data/Wangrong_seismic_data/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
                print 'file not existing' , jJul, Hour 
                continue
            Read_event(Year,jJul,Hour,12, Station, True)
         
    """
    #0. Uniformise file format for file reading ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    if int(Second)>=3600 : 
        print Second
        h = int(Second)//3600
        Hour=int(Hour)+h
        Second = int(Second)%3600
        if int(Hour)>23 :
            jJul = str(jJul+ 1)
            Hour = str(int(Hour-24))
        else : 
            jJul = str(int(jJul))
        if Hour<10 : 
            Hour = '0'+str(Hour)
                
        if int(jJul)<10:
            jJul = '00'+str(jJul)
        elif int(jJul)<100:
            jJul = '0'+str(jJul) 
   
    #1. reading files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st = Stream()
    trE = read('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHE.SAC'%(Year,jJul,Station,Year,jJul,str(Hour)))[0]   
    trE.stats.calib = 1    
    st.append(trE)  
    trN = read('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour)))[0]
    trN.stats.calib = 1    
    st.append(trN)  
    trZ = read('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHZ.SAC'%(Year,jJul,Station,Year,jJul,str(Hour)))[0]
    trZ.stats.calib = 1    
    st.append(trZ)  
    #2. Plotting waveform~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if plot == True : 
        fig = plt.figure()
        fig.canvas.set_window_title('Waveform_%s_%s'%(jJul, Hour))
        st.plot(fig=fig)   
    return st

    
def Stream_Correctionall(st, plot): 
    """return a st whitening and corrected from instrumental response with calibration = 1!!!
    * inputs :
        - st :type str;  stream to be coirrected
        - station :  type str,  station considered
        - plot : if plot == True, plot of un corrected and corrected waveform will be show
    * output : 
        - st : type : stream ; Stream containing the 3 traces ENZ whitening and corrected from instrumental response 
    * exemple : 
      stCorrected = Stream_Correction(st, '1', False)
    """
        
    # 0. calib of traces of stream = 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    l=len(st)
    st_copy = st.copy()
    for i in xrange(l):
        st[i].stats.calib = 1
    st_copy = st.copy()
    stcorrect=Stream()
    #1. Remove mean~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st.detrend('demean') #Mean of data is subtracted 
    #2. Detrend signal with simple ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st.detrend('simple') #Subtracts a linear function defined by first/last sample of the trace
    
    # 3. Correction instrumentale response ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for tr in st : 
        Station = tr.stats.station.split('0')[1]
        if Station in ["1","2","3","4","5","6"]:
           zeros = np.array([0.0+1j*0.0,0.0+1j*0.0,-392.0+1j*0.0,-1960.0+1j*0.0,-1490.0+1j*1740.0,-1490-1j*1740.0])
           poles = np.array([-0.03691+0.03702j,-0.03691-0.03702j,-343.080+1j*0.0,-370.0+467.0*1j,-370.0-467.0*1j,-836.0+1522.0*1j,-836.0-1522.0*1j,-4900.0+4700.0*1j,-4900.0-4700.0*1j,-6900.0+0.0j,-15000+0.0*1j])
           paz_st = {'poles':poles,'zeros': zeros,'gain': 4.3449E17,'sensitivity': 6.1358E+08}
        if Station =='7':
             paz_st = {
            'poles': [-1.48597E-1+1j*1.48597E-1,-1.48597E-1+1j*-1.48597E-1,-2.46936E3+1j*0.0,-4.70635E1+1j*0.0,-3.36765E2+1j*-1.36655E2,-3.36765E2+1j*1.36655E2],
            'zeros': [-3.16174E1+1j*0.0,0.0+1j*0.0,0.0+1j*0.0],
            'gain': 4.8053e+8,
            'sensitivity': 2.4347e+9} 
            # AD facteur = 16.53volt/2^24 sensib=Sv*G(=1)/AD, 
        
        tr.simulate(paz_remove = paz_st,water_level= 1E-4,simulate_sensitivity= False)
        stcorrect.append(tr)
    if plot==True : 
        fig = plt.figure()
        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') 
        plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right='off',      # ticks along the bottom edge are off
        left='off',         # ticks along the top edge are off
        labelleft='off') 
        plt.title('Uncorrected')
        st_copy.plot(fig=fig, label = 'no corrected')
        fig1 = plt.figure()
        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') 
        plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right='off',      # ticks along the bottom edge are off
        left='off',         # ticks along the top edge are off
        labelleft='off') 
        plt.title('Corrected')
        stcorrect.plot(fig=fig1,label = 'corrected')
    return stcorrect
        
    
def Stream_Correction(st, Station, plot): 
    """return a st whitening and corrected from instrumental response with calibration = 1!!!
    * inputs :
        - st :type str;  stream to be coirrected
        - station :  type str,  station considered
        - plot : if plot == True, plot of un corrected and corrected waveform will be show
    * output : 
        - st : type : stream ; Stream containing the 3 traces ENZ whitening and corrected from instrumental response 
    * exemple : 
      stCorrected = Stream_Correction(st, '1', False)
    """
        
    # 0. calib of traces of stream = 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    l=len(st)
    st_copy = st.copy()
    for i in xrange(l):
        st[i].stats.calib = 1
    st_copy = st.copy()
    #1. Remove mean~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st.detrend('demean') #Mean of data is subtracted 
    #2. Detrend signal with simple ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st.detrend('simple') #Subtracts a linear function defined by first/last sample of the trace
    # 3. Correction instrumentale response ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if Station in ["1","2","3","4","5","6"]:
       zeros = np.array([0.0+1j*0.0,0.0+1j*0.0,-392.0+1j*0.0,-1960.0+1j*0.0,-1490.0+1j*1740.0,-1490-1j*1740.0])
       poles = np.array([-0.03691+0.03702j,-0.03691-0.03702j,-343.080+1j*0.0,-370.0+467.0*1j,-370.0-467.0*1j,-836.0+1522.0*1j,-836.0-1522.0*1j,-4900.0+4700.0*1j,-4900.0-4700.0*1j,-6900.0+0.0j,-15000+0.0*1j])
       paz_st = {'poles':poles,'zeros': zeros,'gain': 4.3449E17,'sensitivity': 6.1358E+08}
    if Station =='7':
         paz_st = {
        'poles': [-1.48597E-1+1j*1.48597E-1,-1.48597E-1+1j*-1.48597E-1,-2.46936E3+1j*0.0,-4.70635E1+1j*0.0,-3.36765E2+1j*-1.36655E2,-3.36765E2+1j*1.36655E2],
        'zeros': [-3.16174E1+1j*0.0,0.0+1j*0.0,0.0+1j*0.0],
        'gain': 4.8053e+8,
        'sensitivity': 2.4347e+9} 
        # AD facteur = 16.53volt/2^24 sensib=Sv*G(=1)/AD, 
    
    st.simulate(paz_remove = paz_st,water_level= 1E-4,simulate_sensitivity= False)
    if plot==True : 
        fig = plt.figure()
        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') 
        plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right='off',      # ticks along the bottom edge are off
        left='off',         # ticks along the top edge are off
        labelleft='off') 
        plt.title('Uncorrected')
        st_copy.plot(fig=fig, label = 'no corrected')
        fig1 = plt.figure()
        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') 
        plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right='off',      # ticks along the bottom edge are off
        left='off',         # ticks along the top edge are off
        labelleft='off') 
        plt.title('Corrected')
        st.plot(fig=fig1,label = 'corrected')
    return st
    
    
    
def SelectTracesNEZ(st):
    """returns horizontal traces TrE and trN  and vertical trZ
    * input 
        - st : stream containing 3 components ENZ
    * outputs : 
        - trN, trE, trZ : traces North, east and vertical 
    * exemple: 
        trN, trE, trZ = SelectTracesNEZ(st)
        """
    trN= st.select(component = "N")[0]
    trE= st.select(component = "E")[0]
    trZ= st.select(component = "Z")[0]
    return trN, trE, trZ



def SelectTracesRTZ(st):
    """returns horizontal traces TrR and trT, and vertical TRz 
      * input 
        - st : stream containing 3 components RTZ
    * outputs : 
        - trR, trT, trZ : traces radiale, transverse and vertical 
    * exemple: 
        trR, trT, trZ = SelectTracesRTZ(st)
        """
    trR= st.select(component = "R")[0]
    trT= st.select(component = "T")[0]
    trZ= st.select(component = "Z")[0]
    return trR, trT, trZ
    
    
    
def Rotation(st,angle, plot):
    """return stream rotate thnaks to the Back azimuth between the event and the station
    * input: 
       - st : stream with the 3 components E,N,Z
       - plot : if True, plots of filtered and none filtered will be shown
       - angle : float Bck azimuth between Eevent and Station between [-180,180]
    * output : 
        - stRot : stream with the 3 components R,T,Z
    * exemple : 
       stRot = (st,-179, plot) """
    trN, trE, trZ = SelectTracesNEZ(st)
    trEdata= trE.data
    trNdata= trN.data
       
    if len(trNdata) != len(trEdata):
        raise TypeError("North and East component have different length.")
    if angle < 0:
         angle = 360+angle
    R = trEdata * math.sin((angle + 180.0) * 2 * math.pi / 360.) + trNdata * math.cos((angle + 180.) * 2 * math.pi / 360.)
    T = trEdata * math.cos((angle + 180.) * 2 * math.pi / 360.) - trNdata * math.sin((angle + 180.) * 2 * math.pi / 360.)   
    trR = Trace(R,trN.stats)
    trR.stats.channel = 'BHR'
    trT = Trace(T,trE.stats)
    trT.stats.channel = 'BHT'
    stRot = Stream(traces = [trT,trR,trZ])
    if plot ==True : 
        fig=plt.subplot(211)
        st.plot(fig=fig)
        fig1= plt.subplot(212)
        strot.plot(fig=fig1)
    return stRot



def RMS(st, plot): 
    """return quadradic sum with all signal multiplied by their signs
    * input: 
        - st : stream with the 3 components RTZ or NTZ idem !!
        - plot : if True, plots of filtered and none filtered will be shown
    * output : 
        - tSumtrace : type : trace;  trace signal non decomposed in components
    *exemple
        tSumtrace =  RMS(st, plot)
        ==>  tSumtrace gets all the attributes as a trace"""
        
    trR, trT, trZ = SelectTracesRTZ(st)
    SignR = np.sign(trR.data)
    SignT = np.sign(trT.data)
    SignZ = np.sign(trZ.data)

    SommeQsigne = trR.data**2 + trT.data**2 +trZ.data**2 
   
    SommeSigne = np.sqrt(SommeQsigne) 
    tSumtrace = Trace(SommeSigne,trT.stats)
    tSumtrace.stats.channel = "Sum_quadratic"
    if plot == True :
        fig=plt.figure()
        tSumtrace.plot(fig=fig)
    return tSumtrace
    
   
def RealafterrottrH(st,plot):
    """return quadradic average of the horizontal components 
    * input: 
        - st : stream with the 3 components RTZ !!
        - plot : if True, plots of filtered and none filtered will be shown
    * output : 
        - tSumtrace : type : trace;  trace signal non decomposed in components
    *exemple
        trH =  horizontal quadratic average 
        trZabs = abs average 
        ==> trH gets all the attributes as a trace """
        
    trR,trT,trZ = SelectTracesRTZ(st)
    trZabs=np.abs(trZ.data)
    trZabs=Trace(trZabs,trZ.stats)
    trH  = np.sqrt(trR.data**2 + trT.data**2)/np.sqrt(2)
    trH= Trace(trH,trR.stats)
    trH.stats.channel = "Horizontal"
    if plot == True :
        fig=plt.figure()
        trH.plot(fig=fig)    
    return trH, trZabs
 
    
    
def RealtrH(st,plot):
    """return quadradic sum of the horizontal components with all signal multiplied by their signs
    * input: 
        - st : stream with the 3 components RTZ or NTZ idem !!
        - plot : if True, plots of filtered and none filtered will be shown
    * output : 
        - tSumtrace : type : trace;  trace signal non decomposed in components
    *exemple
        trH =  RMS(st, plot)
        ==> trH gets all the attributes as a trace """
        
    trN,trE,trZ = SelectTracesNEZ(st)
  
    trH  = np.sqrt(trN.data**2 + trE.data**2)/np.sqrt(2)
    trH= Trace(trH,trN.stats)
    trH.stats.channel = "Horizontal"
    trZabs=np.abs(trZ.data)
    trZabs=Trace(trZabs,trZ.stats)
    trZabs.stats.channel = 'Vertical absolute' 
    if plot == True :
        fig=plt.figure()
        trH.plot(fig=fig)    
    return trH,trZabs
        
 
def Stream_PBfilter(st,Fmin, Fmax,plot):
    """ return stream filter with a passband
    * input : 
        - st : streamto filter 
        - Fmin, Fmax  : type, float; min and max frequencies
        - plot :  type, bool; if True, plots of filtered and none filtered are plotted
    * output : 
       -  st_copy : type stream, stream filtered between Fmin and Fmax with a butterworth bandpath
    * exemple : 
     st_filt= Stream_PBfilter(st,0.5, 20,False) """
    #0. copy stream to do not overwrite the raw stream and avoid double filtering~~~~~~~~~~~~~~~~~~~~
    st_copy =st.copy()
    # 1. Pass band filtering  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st_copy.filter('bandpass',freqmin=Fmin,freqmax=Fmax,corners=4,zerophase=True)
    # 2. Plot uncorrected and corrected files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if plot==True : 
        fig = plt.figure()
        plt.title('corrected non filtered_%s-%s'%(Fmin,Fmax))
        st.plot(fig=fig,data_unit="$\\frac{m}{s}$")
        fig1 = plt.figure()
        plt.title('filtered_%s-%s'%(Fmin,Fmax))
        plt.xticks(None)
        st_copy.plot(fig=fig1,data_unit="$\\frac{m}{s}$")
    return st_copy
    
    
    
def LTASTA(st,thres1, thres2,STA, LTA,plotSTA):
    """ return the cut on and off of the LTA/STA list [[cuton, cutoff], [cuton, cutoff]]
    * input : 
        - st : type : stream , stream to filnd STA, LTA
        - thres1 : type; float : cut on limit of STA/LTA values :  after tjis value  the cut on is defined
        - thres2 : type, float : cut off limit of STA/LTA values : after this value the cut off is defined
        - STA : type int :  size of the LTA windows in second : STA = the trace average on this time windows
        - LTA : type int :  size of the LTA windows in second : LTA = the trace average on this time windows
        - plotSTA: type, bool; it true, the trace and it's characteristic function are plotted
        
        RQ: AFTER TESTING IT'S SEEEMS GOOD TO HAVE A RATIO WSTA/WLTA > 1/4 AND A CUT OFF HIGHER THAN CUT ON
    * outputs
        - L_onoff: type np, array : 2D  array of cut on and cut_off time  in number of sample ie time* df where df is the sampling rate [[cuton,cutoff], [cuton1, cutoff1]]
    
    exemple: 
        st = Read_event('15','206','15','1', '1', True)
        stcorrec = Stream_Correction(st,  '1', False)
        stfiltered= Stream_PBfilter(stcorrec,0.5, 20,False)
        LTASTA(stfiltered,2, 2.5,300,1400,True)
    """
    #0. sampling rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    df = st.select(component = "Z")[0].stats.sampling_rate
    #1. characteristic function of the trace following classical LTA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #cft = classicSTALTA(st.select(component = "Z")[0].data, int(STA * df), int(LTA* df))
     #2. characteristic function of the trace following recursive LTA
    cft2 =recursive_sta_lta(st.select(component = "Z")[0].data, int(STA * df), int(LTA* df))
    #3. list of [cuton, cutoff] time in number of samples~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    #max_len = maximum lenght of the triggered event in sample, 
    #max_len_delete = Do not write events longer than max_len into report file.
    L_onoff = trigger_onset(cft2, thres1, thres2, max_len=9e+99, max_len_delete=False)
    #4. plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if plotSTA==True : 
        
        #plotTrigger(st.select(component = "Z")[0], cft, thres1,thres2)
     
        plotTrigger(st.select(component = "Z")[0], cft2, thres1,thres2)
        plt.title('recursive')
    return np.array(L_onoff)
    
def LTASTAtr(tr,thres1, thres2,STA, LTA,plotSTA):
    """ return the cut on and off of the LTA/STA list [[cuton, cutoff], [cuton, cutoff]]
    * input : 
        - tr : type : trace , stream to filnd STA, LTA
        - thres1 : type; float : cut on limit of STA/LTA values :  after tjis value  the cut on is defined
        - thres2 : type, float : cut off limit of STA/LTA values : after this value the cut off is defined
        - STA : type int :  size of the LTA windows in second : STA = the trace average on this time windows
        - LTA : type int :  size of the LTA windows in second : LTA = the trace average on this time windows
        - plotSTA: type, bool; it true, the trace and it's characteristic function are plotted
        
        RQ: AFTER TESTING IT'S SEEEMS GOOD TO HAVE A RATIO WSTA/WLTA > 1/4 AND A CUT OFF HIGHER THAN CUT ON
    * outputs
        - L_onoff: type np, array : 2D  array of cut on and cut_off time  in number of sample ie time* df where df is the sampling rate [[cuton,cutoff], [cuton1, cutoff1]]
    
    exemple: 
        st = Read_event('15','206','15','1', '1', True)
        stcorrec = Stream_Correction(st,  '1', False)
        stfiltered= Stream_PBfilter(stcorrec,0.5, 20,False)
        LTASTA(stfiltered,2, 2.5,300,1400,True)
    """
    #0. sampling rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    df = tr.stats.sampling_rate
    #1. characteristic function of the trace following classical LTA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #cft = classicSTALTA(tr.data, int(STA * df), int(LTA* df))
     #2. characteristic function of the trace following recursive LTA
    cft2 =recursive_sta_lta(tr.data, int(STA * df), int(LTA* df))
    #3. list of [cuton, cutoff] time in number of samples~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    #max_len = maximum lenght of the triggered event in sample, 
    #max_len_delete = Do not write events longer than max_len into report file.
    L_onoff = trigger_onset(cft2, thres1, thres2, max_len=9e+99, max_len_delete=False)
    #4. plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if plotSTA==True : 
        
        ##plotTrigger(tr[0], cft, thres1,thres2)
     
        plotTrigger(tr[0], cft2, thres1,thres2)
        plt.title('recursive')
    return np.array(L_onoff)
    
    
def Trim(st,Second,thres1, thres2,STA, LTA,plotSTA,plotTrim):
    
    """ return the stream cut from the P therical arrival time, to the cut off lTA/ STA are compute on Z component for the moment 
    * inputs : 
        - st : type : stream , stream to find STA, LTA : 
        - Second : type float; theoric arrival time in second of the P wave record in the Event_stations_Caracteristics.txt file
        - thres1 : type; float : cut on limit of STA/LTA values :  after tjis value  the cut on is defined
        - thres2 : type, float : cut off limit of STA/LTA values : after this value the cut off is defined
        - STA : type int :  size of the LTA windows in second : STA = the trace average on this time windows
        - LTA : type int :  size of the LTA windows in second : LTA = the trace average on this time windows
        - plotSTA: type, bool; it true, the trace and its characteristic function are plotted
        - plotTrim:type, bool; it true, the trim stream is plotted
    * ouputs : 
        - st_copy : type, stream; stream trim from the P therical arrival time to the cut off calculated with the recursive LTA/STA methods
   * exemple : 
      st_copy= Trim(stfiltered,'1',60*8,2, 0.5,10,40,True,True)
      ==> 2 plots STA/LTA and the trimed stream
          st trim
      
   """
   #0.sampling rate and time 
    df = st.select(component = "Z")[0].stats.sampling_rate
    start = st.select(component = "Z")[0].stats.starttime
    #1. first trim: keep signal after the P arrival miunus 0.5 second ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st_copy=st.copy()
    st_copy.trim(starttime=start+Second-0.5)
    #2. Filtering of the stream betweenn 0.5 and 15 hz in other to remove a 
    #big part of the noise for the STA/LTA calculation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st_filt = st_copy.copy()
    st_filt= Stream_PBfilter(st,0.5,50,False)
    #3. LTA/STA find time off signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    L_on_off = LTASTA(st_filt,thres1, thres2,STA, LTA,plotSTA)
    off = Second +90 # set by default the off to 90 second
    if len(L_on_off)>0:
        for i in xrange(len(L_on_off)):
           if L_on_off[i][1]/df>Second : 
               #take the first cutoff define
               off= L_on_off[i][1]/df
        
               break
 
    #4. Second signal trim~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    st_copy.trim(endtime = start+off)
    
    if plotTrim==True : 
        fig = plt.figure()
        st_copy.plot(fig=fig)
    return st_copy


def Trimtr(tr,Second,thres1, thres2,STA, LTA,plotSTA,plotTrim):
    
    """ return the stream cut from the P therical arrival time, to the cut off lTA/ STA are compute on Z component for the moment 
    * inputs : 
        - tr : type : trace , trace to find STA, LTA : 
        - Second : type float; theoric arrival time in second of the P wave record in the Event_stations_Caracteristics.txt file
        - thres1 : type; float : cut on limit of STA/LTA values :  after tjis value  the cut on is defined
        - thres2 : type, float : cut off limit of STA/LTA values : after this value the cut off is defined
        - STA : type int :  size of the LTA windows in second : STA = the trace average on this time windows
        - LTA : type int :  size of the LTA windows in second : LTA = the trace average on this time windows
        - plotSTA: type, bool; it true, the trace and its characteristic function are plotted
        - plotTrim:type, bool; it true, the trim stream is plotted
    * ouputs : 
        - st_copy : type, stream; stream trim from the P therical arrival time to the cut off calculated with the recursive LTA/STA methods
   * exemple : 
      tr_copy= Trimtr(stfiltered,'1',60*8,2, 0.5,10,40,True,True)
      ==> 2 plots STA/LTA and the trimed stream
          tr trim
      
   """
   #0.sampling rate and time 
    df = tr.stats.sampling_rate
    start = tr.stats.starttime
    #1. first trim: keep signal after the P arrival miunus 0.5 second ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tr_copy=tr.copy()
    tr_copy.trim(start+Second)
    #2. Filtering of the stream betweenn 0.5 and 15 hz in other to remove a 
    #big part of the noise for the STA/LTA calculation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tr_filt = tr_copy.copy()
    tr_filt= Stream_PBfilter(tr_filt,0.5, 50,False)
    #3. LTA/STA find time off signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    L_on_off = LTASTAtr(tr_filt,thres1, thres2,STA, LTA,plotSTA)
    off = 90+ Second# set by default the off to 90 second
    if len(L_on_off)>0:
        for i in xrange(len(L_on_off)):
           if L_on_off[i][1]/df>int(Second) : 
               off= L_on_off[i][1]/df
               break
        
    #4. Second signal trim~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    tr_copy.trim(endtime = start+off)
    
    if plotTrim==True : 
        fig = plt.figure()
        tr_copy.plot(fig=fig)
    return tr_copy



#st = Stream()
#st = Read_event('15','206','15','1', '1', True)
#stcorrec = Stream_Correction(st,  '1', False)
#stfiltered= Stream_PBfilter(stcorrec,0.5, 20,False)
##LTASTA(stfiltered,2, 2.5,300,1400,True)
#Trim(stfiltered,'1',60*8,2, 0.5,10,40,True,True)
#File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
#B = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
#print ConvertDatestr(B[5])
#Station, Year,jJul, Hour = ConvertDatestr(B[222344])
#st=Read_event(Year,jJul,Hour,12, Station, False)
#Stream_Correction(st, Station, False)
#for j in [1,222344,222345]:
#    
#    Station, Year,jJul, Hour = ConvertDatestr(B[j])
#    Second =B[j][4]
#    print Second
#    print ConvertDatestr(B[j])
#
#   
#    if not  os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
#        print 'file not existing' , jJul, Hour 
#        continue
#    st=Read_event(Year,jJul,Hour,12, Station, True)
#    Trim(st,Second,2, 0.5,10,40,True,True)
