# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 14:25:08 2015

@author: claire
"""

"""Waveform plotting""" 
from Events_caracterisation import *
from Signals_Pre_Processing import *
from obspy.core import  read,UTCDateTime
#from Reading_Seismograms import *
import matplotlib.pyplot as plt
import obspy.signal
import numpy as np
import os.path
import os
from obspy.core.util import gps2DistAzimuth
from obspy.core import read,UTCDateTime,Stream,Trace



"""" List of functions
~~~~~~~~~~~~~~~~~~~~~~
- Plot_Merging_Seismogram :  build a  one day trace merging one hour traces and return the 3 traces 
- Plot_waveform : plot the 3 components of a given event 
- Plot_Merging_2hours_Seismogram: plot 2 hours of signal
- PlotallRot : return radial , transverse and vertical plots for all stations considered for a given event 
"""


#function to merge sesimogram from same stations
def Plot_Merging_Seismograms(Year, jJul):
    """ plot one day (jJul) trace for all the stations pre
    * input : 
        -Year : type str, year of the event ex '15'
        -jJul: type str , jour julian ex '206'
    * output : sent in the jJul file
        - stN, 24h traces North of all the station 1,2,3,4,5,6,7 if files exists
        - stE, 24h traces East of all the station 1,2,3,4,5,6,7 if files exists
        - stZ,24h traces Vertical of all the station 1,2,3,4,5,6,7 if files exists
   ex of use  : 
   ListMtrE,ListMtrN,ListMtrZ =  Plot_Merging_Seismograms(/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,str(1),Year,jJul,Hour))
    #which is on e event in particular
   
   plots all stations all day ! 
   
    """
    Filename = '/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.00.00.00.BHN.SAC'%(Year,jJul,str(1),Year,jJul)
    
    # take the directory name 
    Dirname = os.path.dirname(Filename)
    # Lsit of the file in the directory   
    List_files = os.listdir(Dirname)
    stN = Stream()
    stE= Stream()
    stZ= Stream()
    for i in xrange(0,len(List_files)):

        FileName= List_files[i]
        if 'BHE' in FileName:
            tr = read(Dirname +'/'+FileName)[0]
            tr.stats.calib = 1
            stE.append(tr)
                       
        elif 'BHN' in FileName:
            tr = read(Dirname +'/'+FileName)[0]
            tr.stats.calib = 1
            stN.append(tr)
            
        elif 'BHZ' in FileName : 
            tr = read(Dirname +'/'+FileName)[0]
            tr.stats.calib = 1
            stZ.append(tr)
        
        j = len(List_files)
    #merge all stream with a same component of a given station for all station (station name is given in the metadata)
    stN.merge(method=0, fill_value=None, interpolation_samples=0)
    stE.merge(method=0, fill_value=None, interpolation_samples=0)
    stZ.merge(method=0, fill_value=None, interpolation_samples=0)
    fig1 = plt.figure()
    stN.plot(fig=fig1)
    fig2 = plt.figure()
    stE.plot(fig = fig2)
    fig3 = plt.figure()
    stZ.plot(fig=fig3)
    return stN, stE, stZ

def Plot_Merging_2hours_Seismograms(Year, jJul, hour1, hour2, ListStation):
    """ plot 2 hours (jJul) trace for all the stations pre
    * input : 
        -Year : type str, year of the event ex '15'
        -jJul: type str , jour julian ex '206'
        -hour1 : type str ;first hour 
        -hour2 : type str; second hour
        -ListStation : type list ; list of stations to consider 
    * output : sent in the jJul file
        - stN, 2h traces North of the station in Liststationif files exists
        - stE, 2h traces East of the station in Liststation if files exists
        - stZ,2h traces Vertical of the station in Liststation if files exists
   ex of use  : 
   ListMtrE,ListMtrN,ListMtrZ =  Plot_Merging_Seismograms(/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,str(1),Year,jJul,Hour))
    #which is on e event in particular
       stN, stE,stZ = Plot_Merging_2hours_Seismograms(Year, jJul, '19', '20',  ['1','2','3','4','5','6'])
   plots all stations all day ! 
   
    """

    stN = Stream()
    stE= Stream()
    stZ= Stream()
    for hour in [hour1,hour2]:
        for Station in ListStation:
            if not os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,hour)):    
                print jJul, ' not existing file'              
                continue 
            
            trE = read('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHE.SAC'%(Year,jJul,Station,Year,jJul,hour))[0]
            #trE.stats.calib = 1
            stE.append(trE)
                           
            trN = read('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,hour))[0]
            #trN.stats.calib = 1
            stN.append(trN)
            
            trZ = read('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHZ.SAC'%(Year,jJul,Station,Year,jJul,hour))[0]
           # trZ.stats.calib = 1
            stZ.append(trZ)
     
    #merge all stream with a same component of a given station for all station (station name is given in the metadata)
#    stN.merge(method=-1, fill_value=None, interpolation_samples=0)
#    stE.merge(method=-1, fill_value=None, interpolation_samples=0)
#    stZ.merge(method=-1, fill_value=None, interpolation_samples=0)
    fig1 = plt.figure()
    stN.plot(fig=fig1)
    fig2 = plt.figure()
    stE.plot(fig = fig2)
    fig3 = plt.figure()
    stZ.plot(fig=fig3)
    return stN, stE, stZ       


def Plot_allstation(Composante,event,outfile,Liststation):
    """plot all station of a given event 
    * inputs :
        - Composante : type Str, composante to plot 'N, E, Z
     plotallRot   -event : type list, event = [year, Jul, Hour]
        -outfile : type, str : plot output file
        - Liststion : type list : list of the stations ex ['1','2','3']
    * output :
        - str :  type stream, stream comtaining the traces of the componant and the stations considered 
    * exemple :
       stallN = Plot_allstation('N',[Year, jJul, Hour],False,['1','2','3','4','5','6','7'])
    """
    #ax = plt.subplots(7,sharex=True,sharey=True)
    fig = plt.figure()
    Year = event[0]
    jJul = event[1]
    Hour = event[2]
    L = []
    j=0
    ax11 = plt.subplot(711)
    st =Stream()
    for Station in Liststation :
        if not os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))):    
            print jJul, ' not existing file'              
            continue 
#        
        tr = read('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BH%s.SAC'%(Year,jJul,Station,Year,jJul,str(Hour),Composante))[0]
        tr.stats.calib = 1
        st.append(tr)
       
        print 'Station',  Station
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
    st.plot(fig=fig)
    #if outfile<>False : 
        #st.plot(outfile = outfile, fig=fig)
     
    return st
    
    
def PlotallRot(Year,jJul, Second, Hour, latevent,longevent,depth,LStations):
    """return radial , transverse and vertical plots for all station considered for a given event 
    * input :
        Year :type str  year of the earthquake
        jJul :  type str, julain day of the earthquake 
        Second :  type str ; second of the event 
        Hour: type str; hour of the event consider 
        Latevent: tyope latitu de of the event consider 
        Longevent: type float, longityude of the earthquake condider
        Lstation : type list : list of the station to consider 
    * output : 
        
    * exemple :
    
        """
    b=0
    stR = Stream()
    stT= Stream()
    stZ= Stream()
    
    FileStation = '/home/claire/PHD/Working/Data/station_infos_0.txt'
  
    Catstation=open(FileStation, 'r')
    for LigneSt in Catstation:
     
        if b>0 :
            FLigneSt = LigneSt.split(' ')
           
            ####Station informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      
            Station = FLigneSt[0].split('W0')[1]
            print Station
            if Station in LStations:
                print Station, 'ttt'
                Latst = float(FLigneSt[4])
                Longst = float(FLigneSt[8])
                
                ###Backazimuth, Azzimuth and distance between event and station, radiale distance between sation and event~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                BA, A, Dist, Realdistance, dcartkm =  Azimuth(latevent, longevent, Latst, Longst, depth)
                if not os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))):    
                    print jJul, ' not existing file', Hour, Year, Station            
                    continue 
                
                st=Read_event(Year,jJul,Hour,Second, Station,False)
             
                ###Instrumental correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                stcorrect = Stream_Correction(st, Station, False)
               
                ###Rotation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                strot = Rotation(stcorrect,BA,False)
                Hour = int(Hour)
                ### Trim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                sttrim = strot.trim(starttime = strot[0].stats.starttime + Second, endtime = strot[0].stats.starttime + Second+90)
               
                ###select components~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                trR, trT, trZ=SelectTracesRTZ(sttrim)
                trR.stats.calib = 1
                trT.stats.calib = 1
                trZ.stats.calib = 1
                stR.append(trR)
                stT.append(trT)
                stZ.append(trZ)
        else :
                b+=1 
         
    
    fig1 = plt.figure()
    stR.plot(fig=fig1)
    fig2= plt.figure()
    stT.plot(fig=fig2)
    fig3= plt.figure()
    stZ.plot(fig=fig3)
    return stR, stT, stZ

    
def Plot_oneDay(File):
    plt.figure()
    st = read(File)
    Filename = File.split('/')[len(File.split('/'))-1]
    #tr = Read_sismo(File)[0]
    tr=st[0]
    #tr.plot(type='dayplot')
    t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
    #plt.plot(t,tr.data)
   # plt.show()
    print 'time plot done' 
    
    trFreq= np.fft.fft(tr.data)
    tramplitude = np.abs(trFreq)
    trphase = np.angle(trFreq)
    fs = 100 #sampling rate
    N = len(trFreq)
    xF = trFreq[0:N/2] # since fft is symmetri c over is center 
    fr = np.linspace(0,fs/2,N/2)
    
    """ Amplitude plot """
    plt.subplot(211)
    plt.plot(fr[1:len(fr)-1],np.abs(xF)[1:len(fr)-1])
    plt.ylabel("Amplitude")
    plt.xlabel("Frequency")
    """ Phase plot """ 
    plt.subplot(212)
    plt.plot(fr[1:len(fr)-1],np.angle(xF)[1:len(fr)-1])
    plt.ylabel("Phase")
    plt.xlabel("Frequency")

    """ Inverse fourrier transform : ifft(fft(a))=a frequency domain to time domain""" 
#    plt.subplot(313)
#    trTime = np.fft.ifft(trFreq) 
#    plt.plot(t[0:len(trTime)-1],trTime[0:len(trTime)-1],'r')
#    plt.ylabel("amplitude")
#    plt.xlabel("Time")
    plt.savefig('Plot_Bode/'+Filename.split('B')[0]+'.png')
    return st, fr[1:len(fr)-1],np.abs(xF)[1:len(fr)-1]

def Plot_Filter(tr3, Filename, FilterType,Freqmin,Freqmax,corner):
    
    tr3_filt = tr3.copy()
    tr3_filt.filter('highpass',freq=0.25)
    
    if FilterType=="lowpass" :
        Ft = "LP"
        tr3_filt.filter(FilterType,freq=Freqmin,corners=corner)
        fig = plt.figure()
        tr3_filt.plot(size = (800,600),outfile='Filtered/'+Filename.split('B')[0]+'_Filter_%s_%s.png'%(Ft,str(Freqmin)), fig=fig)
    elif FilterType=='highpass' : 
        Ft = "HP"
        tr3_filt.filter(FilterType,freq=Freqmax,corners=corner)
        fig = plt.figure()
        tr3_filt.plot(size = (800,600),outfile='Filtered/'+Filename.split('B')[0]+'_Filter_%s_%s.png'%(Ft,str(Freqmax)), fig=fig)
    elif FilterType=="bandpass":
        Ft = "BP"
        tr3_filt.filter(FilterType,freqmin=Freqmin,freqmax=Freqmax,corners=corner)
        fig = plt.figure()
        tr3_filt.plot(size = (800,600),outfile='Filtered/'+Filename.split('B')[0]+'_Filter_%s_%s_%s.png'%(Ft,str(Freqmin),str(Freqmax)), fig=fig)
    plt.show()
        
    #tr3.plot(size = (800,600),outfile=Filename.split('B')[0]+'_Filter_%s_%s.png'%(Ft,str(Freq)))
       
    return tr3, tr3_filt


def PlotRecord_session(Files,eq_lat,eq_long):
    """
    """
    #note files should be in a list host is the host site where the data are
    st = stream()
    # Reading the waveforms
    for waveform in Files :
        #add waveform to stream
        st+=read(waveform)
    ## Calculating distance from SAC headers lat/lon
    for tr in st :
        tr.stats.distance = gps2DistAzimuth(tr.stats.sac.stla, tr.stats.sac.stlo,eq_lat, eq_lon)[0]

    #st.filter('bandpass', freqmin=0.1, freqmax=10)
    # Plot
    
    st.plot(type='section', plot_dx=20e3, recordlength=3600, time_down=True, linewidth=.25, grid_linewidth=.25)
    