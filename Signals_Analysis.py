# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 19:23:43 2016

@author: claire
"""
""" Functions List: 

- Envelope : Envelope gives the enveloppe of a signal (trace) and return the sum of this enveloppe
- PeaktoPeak : Peak to Peak : calculate the Peak to peak value of a signal segment and return the time of this peak
- PeaktoPeak_derivative : run peak to peak using the changing sign of the derivative of the signal
- EnvelopeTotal : Total envelope of the signal (high and lower) envelopes calculated using a cubic interpolation of the function between peak and idem between valleys
- HVSR : Horizontal Vertical spectral ratio : if relaively constants there is no site effects
- SmoothSpec : Smooth spectrum of a given trace with the Konno_omachi method
- Saturation : defined if the signal is satured or not
"""


import numpy as np
from obspy.core import read,Trace
from obspy.signal.filter import envelope
from obspy.signal import freqattributes
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing_window,  konno_ohmachi_smoothing, calculate_smoothing_matrix
#from detect_pick import *

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import peakutils
from peakutils.plot import plot as pplot
from scipy.signal import argrelextrema
from Signals_Pre_Processing import *
import copy
def PeaktoPeakPositif(tr):
    """"Peak to Peak : calculate the peak to peak value of a Positive signal value of a signal segment and return the time of this pick """
    #crop the signal on the interval considered     
    tr_copy = tr.copy()
    
    L =[]
    indexesSommet = peakutils.indexes(tr_copy.data,thres = 0.01,min_dist=0)
    indexesVallee = peakutils.indexes(1000-tr_copy.data,thres = 0.01,min_dist=0)

    i=0
    #print 'pt', L 
    if indexesVallee[0]< indexesSommet[0]:
        ptp0 = tr_copy[indexesSommet[0]]-abs(tr_copy[indexesVallee[0]])
        L.append(ptp0)
    for i in xrange(len(indexesSommet)-1):
        Lvalley=[]
        for valley in indexesVallee:
            if valley>indexesSommet[i] and valley<indexesSommet[i+1] :
                Lvalley.append(valley)

            elif valley>indexesSommet[i+1]:
                break
        if len(Lvalley)==1 :

            ptp1 = tr_copy[indexesSommet[i]] - abs(tr_copy[Lvalley[0]])
            L.append(ptp1)
            ptp2 = tr_copy[indexesSommet[i+1]] - abs(tr_copy[Lvalley[0]])
            L.append(ptp2)
        elif len(Lvalley)>1 :
            ptp1 = tr_copy[indexesSommet[i]] - abs(tr_copy[Lvalley[0]])
            L.append(ptp1)
            ptp2 = tr_copy[indexesSommet[i+1]] - abs(tr_copy[Lvalley[len(Lvalley)-1]])
            L.append(ptp2)

    if indexesVallee[len(indexesVallee)-1]>indexesSommet[len(indexesSommet)-1]:
        ptp3 = tr_copy[indexesSommet[i+1]] - abs(tr_copy[valley])
        L.append(ptp3)
    #L = np.array(L)    
    #print len(L),L
    ptpmax = np.amax(L)
    ptpFirst = L[0]

    return ptpmax,ptpFirst


def PeaktoPeak(tr, method, plot):
    """"Peak to Peak : calculate the Pick to pick value of a signal segment and return the time of this pick
     * Input : 
        -tr type :  trace 
        - method :  type int    if 1 : peak to peak is the distance between the peak i and the first valley and the second ptp is the distance between the last valley and the peak i+1
                                if 2 : peak to peak is the distance between the peak i and the minimum valley and the second Peak to peak is the distance between the peak i+1 and the minimum valley. 
        -plot: type boool; if true the plot of peak to peak is plotted 
        
    * output : 
        -MaxptoP:type float; maximum of the peak to peak of the "complet (ie not decompoose in 3 components) of the signal
        -pfirstPtoP type: float; First peak to Peak value : should be first p wave amplitude
    *exemple : 
       ptpmax,ptpFirst =  PeaktoPeak(tr,1,True)
    
    """
    #0. Fuind peaks with a minimum distancebetween peaks in count of 20 and threshold for detcting a peak and valley of 0.001 ~~~~~~~~~~~~~~
    tr_copy = tr.copy()
   
    indexesSommet = peakutils.indexes(tr_copy.data,thres = 0.001,min_dist=20)
    indexesVallee = peakutils.indexes(-tr_copy.data,thres = 0.001,min_dist=20)
    #1. Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if plot==True : 
        time = np.arange(0, tr_copy.stats.npts / tr_copy.stats.sampling_rate, tr_copy.stats.delta)
        plt.figure()
        pplot(time,tr_copy.data,indexesSommet)
        plt.plot(time,tr_copy.data)
        pplot(time,tr_copy.data, indexesVallee)
        plt.show()
   
    i=0 
    #2. Maximum peak to peak calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    L =[] 
    #2.0 Initialisation of the first peak to peak which is equal to the distance betwween the first valley and the first peak. 
    if indexesVallee[0]< indexesSommet[0]:
        ptp0 = tr_copy[indexesSommet[0]]+abs(tr_copy[indexesVallee[0]])
        L.append(ptp0)
    # 3.0 measure peak to peak between 2 peaks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i in xrange(len(indexesSommet)-1):
        Lvalley=[]
        #3.1 record all the valleys indexes that are between peaks i and i+1 to then take into account only the fist and the second peaks into account (maybe better to take the minimum! )
        for valley in indexesVallee:
            if valley>indexesSommet[i] and valley<indexesSommet[i+1] : 
                Lvalley.append(valley)
    
            elif valley>indexesSommet[i+1]:
                break
        #3.2 calculate pêak to peak ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if len(Lvalley)==1 : 
            
            ptp1 = tr_copy[indexesSommet[i]] + abs(tr_copy[Lvalley[0]])
            L.append(ptp1)
            ptp2 = tr_copy[indexesSommet[i+1]] + abs(tr_copy[Lvalley[0]])
            L.append(ptp2)
        elif len(Lvalley)>1 :
            # peak to peak is the distance between the peak i and the first valley and the second ptp is the distance between the last valley and the peak i+1
            if method == 1 : 
                
                ptp1 = tr_copy[indexesSommet[i]] + abs(tr_copy[Lvalley[0]])
                L.append(ptp1)
                ptp2 = tr_copy[indexesSommet[i+1]] + abs(tr_copy[Lvalley[len(Lvalley)-1]])
                L.append(ptp2)
            #peak to peak is the distance between the peak i and the minimum valley and the second Peak to peak is the distance between the peak i+1 and the minimum valley. 
            elif method == 2 :
                Lvaluevalley = []
                if j in Lvalley : 
                    Lvaluevalley.append(abs(tr_copy[Lvalley[j]]))
                Valleymax= max(Lvaluevalley)
                ptp1 = tr_copy[indexesSommet[i]] + max(Lvaluevalley)
                L.append(ptp1)
                ptp2 = tr_copy[indexesSommet[i+1]] + max(Lvaluevalley)
                L.append(ptp2)
                    
                
                
    #4.0 end if the signal finished with a valley the last peak to peak is the distance between the last peak and the last valley ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
    if indexesVallee[len(indexesVallee)-1]>indexesSommet[len(indexesSommet)-1]:
        ptp3 = tr_copy[indexesSommet[i+1]] + abs(tr_copy[valley])
        L.append(ptp3)   
    #5. Take the maximum of the peak to peak ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ptpmax = np.amax(L)
    # 6. the first vallue of the peak to peak which is suppose to be the P wave amplitude ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ptpFirst = L[0]

    return ptpmax,ptpFirst

def Peakindex_derivative(tr,thresMax):
    """ ruturn the index of the peaks that are 0.98 greater than the max peak
    *input : 
        - tr: type : trace , trace of the signal to find max and min 
    *output :
        - LPtP :type list;  Liste of Peak to Peak values
        - MaxPtP :  Maximum of PtP values 
    *exemple :
       LPtP, MaxPtP=  PeaktoPeak_derivative(tr,True)
    """
    #0. Sampling rate
    df = tr.stats.sampling_rate
    #1. Derivative of the signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Dsignal= []
    for i in xrange(len(tr.data)-1) : 
        f1 = (tr.data[i+1]-tr.data[i])/df 
        Dsignal.append(f1)
    #2. index where f1 == 0 ie index of extremumns~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    asign = np.sign(Dsignal)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)

    Zeros = np.argwhere(signchange==1)
    indexmax= np.argmax(tr.data)
    Peakselect = []
    for index in Zeros : 
        if tr.data[index[0]]>=thresMax * tr.data[indexmax]:
      
            Peakselect.append(index[0])
        
    return Peakselect
    
def PeaktoPeak_derivative(tr,plot):
    """ run peak to peak using the changing sign of the derivative of the signal
    *input : 
        - tr: type : trace , trace of the signal to find max and min 
    *output :
        - LPtP :type list;  Liste of Peak to Peak values
        - MaxPtP :  Maximum of PtP values 
    *exemple :
       LPtP, MaxPtP=  PeaktoPeak_derivative(tr,True)
    """
    #0. Sampling rate
    df = tr.stats.sampling_rate
    #1. Derivative of the signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Dsignal= []
    for i in xrange(len(tr.data)-1) : 
        print 'derivative ', len(tr.data)
        f1 = (tr.data[i+1]-tr.data[i])/df 
        Dsignal.append(f1)
    #2. index where f1 == 0 ie index of extremumns~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    asign = np.sign(Dsignal)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    print type(signchange)
    Zeros = np.argwhere(signchange==1)
    print Zeros, 'zeros'
    LPtP=[]
    #3. Value of extremums ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i in xrange(len(Zeros)-1):
        PtP = abs(tr.data[Zeros[i]])+abs(tr.data[Zeros[i+1]])
        LPtP.append(PtP)
    #4. Maximum~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MaxPtP = max(LPtP)
    #plot scatter, trace derivative and trace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if plot ==True : 
        time = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
        plt.figure()
        plt.scatter(time[0:len(signchange)],signchange*100000, marker = '+', color = 'r')
        plt.plot(time, tr.data)
        plt.plot(time[0:len(Dsignal)], Dsignal, 'r')
    return LPtP, MaxPtP
    
#    
def Envelope(tr,plot):
    
    """Envelope gives the enveloppe of a signal (trace) and return the sum of this enveloppe
    * input : 
        - tr : type trace; trace 
        -plot type bool; if True waves and envelope are plotted 
    *output : 
        - Envelope : envelope of the signal 
        - t: type list; time where data are computed
        - SumEnvelope : type float;  Sum of the data values of the enveloppe
        - MaxEnv  ; type float; maximum of the envelope
    *exemple :
        Envelope, t, SumEnvelope = Envelope(tr)
    """
    Envelope = envelope(tr.data)
    npts = tr.stats.npts
    fs = tr.stats.sampling_rate
    t = np.arange(0,npts/fs,1/fs)
    #integration of the envelope
    SumEnvelope = sum(Envelope)
    #max~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EnvMax = max(Envelope)
    indexmax= np.argwhere (Envelope==EnvMax)
    if plot==True: 
        plt.figure()
        plt.plot(t,tr.data )
        plt.plot(t, Envelope, '.g')
    print  EnvMax,indexmax
    return Envelope, t, SumEnvelope, EnvMax
    
def EnvelopeTotal(tr, plot) : 
    """Total envelope of the signal (high and lower) envelopes calculated using a cubic interpolation of the function between peak and idem between valleys
    * input : 
        - tr : type trace; trace 
        -p lot type bool; if True waves and enevlopes are plotted 
    *output :  , q_l, MaxEnv, indexMax
        - q_u : type, array; upper envelope
        - q_l : type array; lower envelope
        - MaxEnv: type float; maximum of the envelope amplitude (max sul upper and abs (lower) envelopes)
        - indexmax : type, int, index of the maximum of the envelope

    *exemple :
        Env_u, Env_l, MaxEnv, indexMax = EnvelopeTotal(tr, True)
    """
   #0. initialisation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    q_u = np.zeros(tr.data.shape) 
    q_l = np.zeros(tr.data.shape) 
    
    # 1. Prepend the first value of (s) to the interpolating values. This forces the model to use the same starting point for both the upper and lower envelope models.
    #upper envelope
    u_x = [0,]
    u_y = [tr.data[0],]
    #lower envelope
    l_x = [0,]
    l_y = [tr.data[0],]
    
    # 2. find peaks and valleys  and mark their location an dvalues in u_x,u_y,l_x,l_y respectively.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for k in xrange(1,len(tr.data)-1):
        if (np.sign(tr.data[k]-tr.data[k-1])==1) and (np.sign(tr.data[k]-tr.data[k+1])==1):
            u_x.append(k)
            u_y.append(tr.data[k])
    
        if (np.sign(tr.data[k]-tr.data[k-1])==-1) and ((np.sign(tr.data[k]-tr.data[k+1]))==-1):
            l_x.append(k)
            l_y.append(tr.data[k])
    #3. Append the last value of (s) to the interpolating values. 
    #This forces the model to use the same ending point for both the upper and lower envelope models.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    u_x.append(len(tr.data)-1)
    u_y.append(tr.data[-1])
    
    l_x.append(len(tr.data)-1)
    l_y.append(tr.data[-1])
    
    #4. Fit suitable models to the data. Here I am using cubic splines, 
    #   similarly to the MATLAB example given in the question.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    u_p = interp1d(u_x,u_y, kind = 'cubic',bounds_error = False, fill_value=0.0)
    l_p = interp1d(l_x,l_y,kind = 'cubic',bounds_error = False, fill_value=0.0)
    
    
    #5. Evaluate each model over the domain of (s)
    for k in xrange(0,len(tr.data)):
        q_u[k] = u_p(k)
        q_l[k] = l_p(k)
     
    #6. Calculate max envelope~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sum = q_u+abs(q_l) 
    MaxEnv= max(Sum)
    indexMax= np.argwhere(Sum == MaxEnv)
    if plot==True: 
        #time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        npts = tr.stats.npts
        fs = tr.stats.sampling_rate
        t = np.arange(0,npts/fs,1/fs)
        #plot
        plt.figure()
        #trace plot
        plt.plot(t, tr.data, 'k')
        # upper and lower envelope plots
        plt.plot(t, q_u,'.g',markersize=1 )
        plt.plot(t, q_l, '.g',markersize=1 )
        #max enveloppe locatiojn and value
        plt.scatter(indexMax/fs, MaxEnv, marker= '+', color= 'g')
    print MaxEnv, indexMax
    return q_u, q_l, MaxEnv, indexMax



def HVSR(SH,SV):
    """ Horizontal Vertical spectral ratio : if relaively constants there is no site effects
    * input : 
        - SH : type array, horizontal spectrum
        - SV: type array , vertical spectrum
    * output : 
        - HVSR :  type array; horizontal/vertical spectrum ratio
    * exemple :
    """
    HVSR = SH/SV
    return HVSR
    
def KonnoSmoothSpec(tr, bandwidth, plot):
    """ Smooth spectrum of a given trace with the Konno_omachi method
    * inputs : 
        - tr: type: trace; seimic trace
        - bandwidth : type int; number of point over wich the spectrum is smoothed, after testing 110 seems good
        - plot : type Bool; if True smooth, and un smooth spectrum are plotted
    * outputs : 
        - fr :  type array; frequencies of the spectrum
        - SmoothSpectre : type array; Smooth spectrum
    * exemple :
    """
    #  0. Spectre of the trace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
    trFreq= np.fft.fft(tr.data, n= 10000)
    fs = tr.stats.sampling_rate
    N = len(trFreq)
    #half spectre since fft is symmetric over is center 
    xF = trFreq[0:N/2] # 
    #frequencies
    fr = np.linspace(0,fs/2,N/2)
    
    #1. smoothing Spectrum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #CAREFULL SMOTTHING GONNA CHANGE A LOT THE SPECTRUM BUT DOING THE SPECTRUML RATIO IS OK!
    SmoothSpectre = konno_ohmachi_smoothing(np.abs(xF), fr, bandwidth=bandwidth, count=1, enforce_no_matrix=False, normalize=False)
    
    if plot== True : 
        plt.figure()
        plt.plot(fr,np.abs(xF), 'r',label = 'non Smooth')
        plt.plot(fr,SmoothSpectre , 'k', label='Smooth')
        plt.ylabel("Amplitude")
        plt.xlabel("Frequency Hz")
        plt.legend()
    return fr, SmoothSpectre    

def Smooth(x,fr,plot,window_len=110, window='flat'):
    """smooth the data using a window with requested size.
        This method is based on the convolution of a scaled window with the signal
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
    
  input:
      - x: the input spectrum
      - window_len: the dimension of the smoothing window; should be an odd integer
      - window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
      -plot :  type Bool; if plot == True, plot of the smooth spectra will be return 
 
    output:
        the smoothed signal
     example:

      x=sin(t)+randn(len(t))*0.1
      y=smooth(x)
    
    see also: 
   
 
      NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
      """ 

    if x.ndim != 1:

        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
         return x     
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
         raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'" 
    #build an array quickly : numpy_r
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    f = np.r_[fr[window_len-1:0:-1],fr,fr[-1:-window_len:-1]]
 
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    f2=np.convolve(w/w.sum(),f,mode='valid')
    if plot== True :
        plt.figure()
        plt.plot(f2, y,'.g')
    return y ,f2

def Spectrum(trf, plot):
    """Amplitude spectrum of a given trace 
    * inputs : 
        - tr: type: trace; seimic trace
        - plot : type Bool; if True smooth, and un smooth spectrum are plotted
        
    * outputs : 
        - fr :  type array; frequencies of the spectrum
        -Spectre : type array; amplitude spectrum
    * exemple :
        fr, spectrum = Spectrum(tr, True)
    """
    #  0. Spectre of the trace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    t = np.arange(0, trf.stats.npts / trf.stats.sampling_rate, trf.stats.delta)
    trFreq= np.fft.fft(trf.data, n=1000) 
    fs = trf.stats.sampling_rate
    N = len(trFreq)
    #half spectre since fft is symmetric over is center 
    xF = trFreq[0:N/2] # 
    #frequencies
    fr = np.linspace(0,fs/2,N/2)
    if plot==True: 
        plt.figure()
        plt.plot(fr,np.abs(xF), 'k')
        plt.ylabel("Amplitude")
        plt.xlabel("Frequency Hz")
        plt.show()
    return fr, np.array(abs(xF))

def SmoothSpectrum(tr, plot):
    """Amplitude spectrum of a given trace 
    * inputs : 
        - tr: type: trace; seimic trace
        - plot : type Bool; if True smooth, and un smooth spectrum are plotted
    * outputs : 
        - fr :  type array; frequencies of the spectrum
        -Spectre : type array; amplitude spectrum
    * exemple :
        fr, spectrum = Spectrum(tr, True)
    """
    #  0. Spectre of the trace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
    #set the niumber of poitn or else each spectrum for a given event won't have the same frquency content and erros will happened
    trFreq= np.fft.fft(tr.data, n=100000)
    fs = tr.stats.sampling_rate
    N = len(trFreq)

    #half spectre since fft is symmetric over is center 
    xF = trFreq[0:N/2] # 

    #frequencies
    fr = np.linspace(0,fs/2,N/2)
    # abs = module d'un nombre complexe sous python!!!!!
    xfsmooth,frsmooth = Smooth(np.array(abs(xF)),fr,plot,window_len=300, window='flat')
    if plot==True: 
        plt.figure()
        plt.plot(fr,np.abs(xF), 'k')
        plt.ylabel("Amplitude")
        plt.xlabel("Frequency Hz")
        
    return frsmooth, np.array(xfsmooth)
    
def Saturationold(tr,plot):
    """ return true and the component if the trace is satured find if there is a plateau formed
    input : 
    tr: type trace , seimic trace to analyse
    output : 
    Stat type bool, if true the signal is considerd as satured"""
    indexesSommet = peakutils.indexes(tr.data,thres = 0.01,min_dist=20)
    indexMax= np.argmax(tr.data)
    
    indexesSommet=np.delete(indexesSommet,np.argwhere((indexesSommet==indexMax)))
    
    if len(indexesSommet)<>0:
        MAx2=np.argwhere(tr.data==np.max([tr.data[i] for i in indexesSommet]))
    else :
        MAx2 =[indexMax]
    if len(MAx2)==0:
        MAx2 =indexMax
    else : 
        MAx2 = MAx2[0]
    
    if indexMax<len(tr.data)-2 and MAx2<len(tr.data)-2:
       
        if tr.data[indexMax+1]==tr.data[indexMax] or tr.data[MAx2+1]==tr.data[MAx2] :
            Sat = True
        else : 
            Sat= False
    else : 
        Sat = True
    return Sat

def Saturation(tr,plot):
    """ return true and the component if the trace is satured find if there is a plateau formed
    input : 
    tr: type trace , seimic trace to analyse
    output : 
    Stat type bool, if true the signal is considerd as satured min_dist min poiunts between 2 peaks"""
    #selection of peaks that have asignal greater than 0.7*np.max(tr.data)
    indexselected = Peakindex_derivative(tr,0.7) 

    indexselected.sort()
    indexselectedbis = copy.copy(indexselected)
    
    #initialisation of the list of saturated peaks index that are greater than 0.95 than tr.data+1
    Lsat =[]
    for i in xrange(len(indexselected)-1) : 
        if tr.data[indexselected[i+1]] == tr.data[indexselected[i]]:
            Lsat.append(True)
        elif tr.data[indexselected[i+1]]>=0.997 * tr.data[indexselected[i]] and tr.data[indexselected[i+1]]<=(1+0.013)*tr.data[indexselected[i]] and abs(indexselected[i+1]-indexselected[i])<3:
            Lsat.append(True)
        else:
            Lsat.append(False)
            indexselectedbis.remove(indexselected[i])
    #Thje siganl is satured if one peak have a neigbourg with an amplitude greater than 0.98 of his amplitude
    a = len(np.where(Lsat)[0])

    if a>=2:
        Sat = True
    else :
        Sat = False
    #plot
    if plot==True : 
        time = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
        plt.figure()
        
        pplot(time,tr.data,indexselected)

        plt.plot(time,tr.data)
        #pplot(time,tr.data,Listindex60)
  
        plt.show()
    return Sat,Lsat

def SNR(tr, trTrim, LimSNR, Second,plot):
    MeanPSignal = np.mean(trTrim.data**2)
    Maxsignal = np.max(trTrim.data)**2
    tr1 = tr.copy()
    trnoise = tr1.trim(starttime =tr1.stats.starttime+Second-190,endtime=tr1.stats.starttime+Second-100)
    MeanPnoise = np.mean((trnoise.data)**2)
    SNR = Maxsignal/MeanPnoise
    #SNR = MeanPSignal/MeanPnoise
    if SNR>LimSNR :
        valid = True
    else :
        valid = False
        
    if plot==True:
        figurenoise =plt.figure()
        
        trnoise.plot(fig=figurenoise,color='b',linewidth=2)
        trTrim.plot(fig=figurenoise,color='r',linewidth=2)
        tr.plot(fig=figurenoise,starttime =tr.stats.starttime ,endtime =tr.stats.starttime +Second +90 )
    return valid, SNR
    
    
def SNRstd(tr, trTrim, LimSNR, Second,plot):
    STDSignal = np.std(trTrim.data)
    tr1 = tr.copy()
    trnoise = tr1.trim(starttime =tr1.stats.starttime+Second-190,endtime=tr1.stats.starttime+Second-100)
    STDnoise = np.std(trnoise.data)
    SNR = STDSignal/STDnoise
    #SNR = MeanPSignal/MeanPnoise
    if SNR>LimSNR :
        valid = True
    else :
        valid = False
        
    if plot==True:
        figurenoise =plt.figure()
        
        trnoise.plot(fig=figurenoise,color='b',linewidth=2)
        trTrim.plot(fig=figurenoise,color='r',linewidth=2)
        tr.plot(fig=figurenoise,starttime =tr.stats.starttime ,endtime =tr.stats.starttime +Second +90 )
    return valid, SNR
#test 
    
#st = read('/home/claire/PHD/Working/Data/Wanrong_array_2015/2015/R206.02/GSW03.15.206.10.00.00.BHN.SAC')
#fig = plt.figure()
#st=Stream_Correction(st, '3', False)
#st.trim(st[0].stats.starttime + 25*60+40, st[0].stats.starttime + 26*60+40)
#PeaktoPeak_derivative(st[0],True)
#Envelope(st[0], True)
#EnvelopeTotal(st[0], True)
#SmoothSpec(st[0], 150, True)
