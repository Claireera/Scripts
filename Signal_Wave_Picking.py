# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 18:23:30 2016

@author: claire
"""

"""Wave pickingphase S and P following (Baillard 2014) methods

Kurtosis :calculate the curtosis of a trace
Cov : calculmate the covariance of the signal
polarisation: calculate the degreee of rectinearility and the dip of maximum polarization
"""
import math
import numpy as np
from Signals_Pre_Processing import *
import scipy
from datetime import timedelta
from heapq import merge
def CFkurtosis(tr,T):
    """return the caracteristic function based on the cumulative kurtosis calculate over a 5 second time window
    input : 
        - tr: trace 
    output :
        - Cf Characteristic function which is the cumulative kurtosis calculated over T second time windows 
    exemple :
        for a period T = 5 seconds 
        Cf = CFkurtosis(tr.data,5) """
        
    #100 sample per seconds thus 5 s =500 samples first one  is 0+500. 
    N = T*100
    Cf =np.zeros((1,N)) 
    cf = []
    for i in xrange(N, len(tr.data)):
        data = tr.data[i-T*100+1:i]
        #Kurtosis calculation as follow the Pearson definition
        Kurtosis =  scipy.stats.kurtosis(data, axis = 0, bias=False, fisher=False)
        cf.append(Kurtosis)
    Cf = np.append(Cf, np.array(cf))
    #Cf = np.append(np.zeros(N),Cf)
    return Cf
    
def F2(Cf):
    """ remove the negative slope of the characteristic function calculated as the sliding kurtosis fonction
    input :
        Cf: 1D array, Kurtosis characteristic funtion calculated with CFkurtosis 
    output:
        F2: 1D array, removal of negative slopes of Cf
    exemple :
        F2 = F2(Cf)
    """
    F2=np.zeros(Cf.shape)
    F2[0]=Cf[0]
    print Cf.shape, "cf"
    for i in xrange(len(Cf)-1):
        df1 = Cf[i+1]-Cf[i]
        if df1<0:
            delta =0
        else :
            delta =1
            
        F2[i+1] =F2[i]+delta*df1
    return F2
        
def F3(F2):
    """remove the lienar trendfrom F2
    input:
        - F2: 1D array, removal of negative slopes of kurtosis characteristic function of a given trace
    output:
        - F3: 1D array, removal of the linear trend of F2
    exemple:
        Cf = CFkurtosis(tr.data,5)
        F2 = F2(Cf)
        F3 = F3(F2)
        """
    
    F3= np.zeros(F2.shape)
    a = (F2[len(F2)-1]-F2[0])/(len(F2)-1)
    b = F2[1]
    for i in xrange(1,len(F3)):
        F3[i]=F2[i]-(a*(i-1)+b)
    #   shift the F3 component to the left
    F3=np.array([F3[i-1] for i in xrange(1,len(F3))])
    F3 = np.append([0],F3)
    return F3
    
def MeanF3(F3,WT):
    """
    imput:
        -F3: 1D array, removal of the linear trend of F2
        -WT: int, smoothing windowss size in second 
    output :
        -F3prim: 1D array, F3 smoothed
        
    exemple:
        F3prim = MeanF3(F3,1)
"""    
    F3prim=np.zeros(F3.shape)
    Wt=int(100*WT)
    for i in xrange(len(F3)-Wt):
        F3prim[i]=np.mean(F3[i:i+Wt])
        
    return F3prim
    
def F4(F3):
    """final transformation:
    input : 
        -F3: 1D array,
     
    output :
        -F4: 1D array,picking minima ion F4 ll give a good estimate on phase onset 
    exemple :
    F4 = F4(F3)
    """
    F4 = np.zeros(F3.shape)
    T=np.zeros(F3.shape)
    for i in xrange((len(F3)-1)):
        T[i]=F3[i]-max(F3[i],F3[i+1])
        if T[i]<0:
            F4[i]=T[i]
        else : 
            F4[i]=0
    return F4
    
def F4prim(F4):
    """final transformation:
    input : 
        -F3: 1D array,
     
    output :
        -F4: 1D array,picking minima ion F4 ll give a good estimate on phase onset 
    exemple :
    F4 = F4(F3)
    """
    F4prim = np.zeros(F4.shape)
    T=np.zeros(F4.shape)
    for i in xrange((len(F4prim)-1)):
        # f4 IS NEGATIVE
        F4prim[i]=F4[i]-max(F4[i],F4[i+1])
        
    return F4prim
    
def Covariance(st):
    """ calculate the covariance of a signal to determine the polarisation
    input: 
        - st, stream containing the 3 trace (E,N, Z)
    output: 
        - C, covariance matrix
    exemple"""
    trN,trE,trZ= SelectTracesNEZ(st)
  
    #Covariance matrix
    size= len(trN.data)
    c11= np.dot(trN.data, trN.data)/size
    c12= np.dot(trN.data, trE.data)/size
    c13= np.dot(trN.data, trZ.data)/size
    c21 = np.dot(trE.data, trN.data)/size
    c22 = np.dot(trE.data, trE.data)/size
    c23 = np.dot(trE.data, trZ.data)/size
    c31 = np.dot(trZ.data, trN.data)/size
    c32 = np.dot(trZ.data, trE.data)/size
    c33 = np.dot(trZ.data, trZ.data)/size
    Cov = np.array([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])
    return Cov
    
def eignevalues(st):
    """determines the eiugenvalues of the sigal
    inputs: 
        - trN,trE,trZ,1D arrays data trace N, E and Z
    outputs:
        -w: order eigenvalues
        -v: order eigenvectors V[:,i] is the eigen vector corresponding to the eignevalues w[i]
    exemple"""
    #trN,trE,trZ= SelectTracesNEZ(st)
    # concatenate the matrix: column concatenations c1 = trN , c2= trE, c3= trZ
    Cov= Covariance(st)
    w, v = np.linalg.eig(Cov)
    if w<>[]:
        # argument of the max of eigen *
        i1 = np.argmax(w)
        i3 = np.argmin(w)
        arg= [0,1,2]
        arg.remove(i1) 
        arg.remove(i3)
        i2 = arg[0]
        w= np.array([w[i1],w[i2],w[i3]])
        v = np.array([v[:,i1],v[:,i2],v[:,i3]])
    
    return w, v
    
def Polarisation(L1,L2,L3,U1):
    """ calculate the degreee of rectinearility and the dip of maximum polarization
    L and U are defined C.Ui =LiUi 
    input: 
    L1: float, first eigenvalue as L1>L2>L3
    L2: float, Second eigenvalue as L1>L2>L3
    L3 :float, Third eigenvalue as L1>L2>L3
    U1:1D array,unit first eigenvector (ax of the maximum polarisation) 
    output:
    Rectinearility : float, degree of rectinearility define following (Jurkevics 1988) Rec = 1: rectilinear polarisation, Rec =Circular polarisation    
    Dip : float, the dip of maximum polarisation (Vidale 1986) in degrees
    exemple:
    """
    Rectinearility =1-((L2+L3)/(2*L1))
    Dip =math.degrees(math.atan(U1[2]/math.sqrt(U1[1]**2+U1[0]**2))) #result given in degrees
    return Dip, Rectinearility
    
def P_S_Characterisation(st,Twin,alpha, k, Tp,startsecond, endsecond): 
    """ diffrenciate surface and body waves and P and S onset (P is suupose to have mostly vertical motion and S wave horizontakl one (incident waves have small incidence angle ))
    polarisation is calculated on 10s moving windows with 20% overlaps following the Vidale 1986 method between 0.02-15 Hz"     
    input :
        - st : stream containing TrN, TRE and TrZ
        - Twin : float, 4s time of the sliding windows see baillard but should cgahnge 
        - Alpha: float, parameter in dip rectilinearity
        - k: int, indice of the coordinate of the pic identify with F4
    output :
        onset; string, nature of the pic identify P, S or None which is not a body wave   
        Dr: 1D array , dip reactilinearity function
    exemple :"""
    
    #1. dip-rectinearility fonction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    win_len = 1 # 10 seconds of sliding windows lenght following engineering speedy technique to evaluate Seiemsic site effects
    win_frac = 0.2 #fraction of sliding windows step following idem 
    frqlow = 0.02 #lower frq for pm
    frqhigh = 15 #higher frq fo pm

    # dip and rectinearility
    i= 0
    #initailisation of the diprectilinearity function
    DR=[]
    trN,trE,trZ= SelectTracesNEZ(st)
    #number of sample per second
    Nsample = trN.stats.sampling_rate
    J= int(Twin*Nsample)
    
    #calculation of rectinearility and dip over the time on a sliding windows of size Twin 
    for i  in xrange(len(trZ.data)-J):
        stcopy = st.copy()
        stcopy.trim(starttime=trN.stats.starttime+ i/Nsample, endtime=trN.stats.starttime + (i+J)/Nsample)     
            
        
        if len(trN.data)==J+1:
            trN,trE,trZ= SelectTracesNEZ(stcopy)        
            print trN.stats.starttime, "f4 strat time fr" 
            v,w = eignevalues(stcopy)
            if v<>[]:
                L1 = v[0]
                L2 = v[1]
                L3 = v[2]
                U1 = w[0] 
                #rectilinearity and Dip
                dip,rect = Polarisation(L1,L2,L3,U1)
                #dip rectinearility at i
                dr =  rect * (np.sign(alpha*math.sin(dip)-rect))
                DR.append(dr)
            else :
                DR.append(0)
        else :
            DR.append(0)
    print len(DR), "dr-f4"
    if DR[k]>Tp :
        onset ='P'
    elif DR[k]<-Tp :
        onset ='S'
    else : 
        onset = 'None'
    return DR, onset 


def WavePicking(tr,T, Second,distance,plot):
    """  give start of p and s 
    -input:
        -tr : trace 
        -T : float, 5s time of the sliding windows see baillard but could be change trial shoaw goiod result with 5s. 
        plot: bool, if true plot will be draw
    -output
        - minF4p: float coordiantes of the p start in second
        - minF4s: float coordiantes of the s start in second
    -exemple:
        minF4p, minF4s = WavePicking(tr,5,False)"""
    tr_copy = tr.copy()
    #NESURTOUT PAS FILTERER
    #tr_copy.filter('bandpass', freqmin=0.2, freqmax= 15)
    if Second > 60*5: 
        tr_copy.trim(starttime =tr.stats.starttime + Second-60*5)
#    if 3600 > Second+60*5:
#        tr_copy.trim(endtime =tr.stats.starttime + Second+60*10)
    if tr.stats.starttime +SecondP + P_Sdelay+120<tr.stats.endtime:
        tr_copy.trim(endtime =tr.stats.endtime + SecondP+P_Sdelay+120)

    Cf = CFkurtosis(tr_copy,T) 

    f2 = F2(Cf)

    f3 = F3(f2)

    f3prim = MeanF3(f3,0.5) #sliding windows of 0.5 s show good results 
 
    f4 = F4(f3prim)
    f4prim = F4prim(f4)
    # P and S coordinates 
    minF4p = np.argmin(f4)
    print len(f4),'lenf4', minF4p, 'p arrival'
    #minimum of points between p and s start function of distance  (calculated with vs =2.6 km.s-1 and vp = 4.1km.s-1)
    if distance<=20 :
        nbPtmin= 20
    
    elif 20<distance<=40 :
        nbPtmin= 50
    elif 40<distance<=60 : 
        nbPtmin= 60
    elif 60<distance<=80 :
        nbPtmin= 90
    elif 80<distance<=100 : 
        nbPtmin= 110
    else:
        nbPtmin= 140

    if minF4p+3000<len(f4):
        min2F4p = np.argmin(f4prim[minF4p+nbPtmin:minF4p+3000])
        minF4s = minF4p+nbPtmin+min2F4p# strat 100 point after the first minimum found to find the second one 
    else:
        print 'CLOSE TO THE END OF THE SIGNAL'
        if minF4p+nbPtmin<len(f4):
            min2F4p = np.argmin(f4prim[minF4p+nbPtmin:len(f4)])
            minF4s = minF4p+100+min2F4p
        else:
            min2F4p = np.argmin(f4prim[minF4p+1:len(f4)])
            minF4s = minF4p+1+min2F4p
            'SIGNIAL TOO SHORT'
    if plot == True:
        component = tr.stats.channel
        f, (ax0,ax1, ax2, ax3, ax4,ax5) = plt.subplots(6, sharex=True)
        x = range(len(tr_copy.data))
        ax0.plot(x,tr_copy.data)
        ax0.axvline(x=minF4p,color= 'r')
        ax0.axvline(x=minF4s,color= 'aqua')
        ax0.set_ylabel('$velocity m.s^-1$')
        ax1.plot(x,Cf)
        ax1.axvline(x=minF4p,color= 'r')
        ax1.axvline(x=minF4s,color= 'aqua')
        ax1.set_ylabel('$Cf$')
        ax2.plot(x, f2)
        ax2.axvline(x=minF4p,color= 'r')
        ax2.axvline(x=minF4s,color= 'aqua')
        ax2.set_ylabel('$F2$')
        ax3.plot(x,f3prim)
        ax3.set_ylabel('$F3$')
        ax3.axvline(x=minF4p,color= 'r')
        ax3.axvline(x=minF4s,color= 'aqua')
        ax4.plot(x, f4)
        ax4.axvline(x=minF4p,color= 'r')
        ax4.axvline(x=minF4s,color= 'aqua')
        ax4.set_ylabel('$F4$')
        ax5.plot(x, f4prim)
        ax5.axvline(x=minF4p,color= 'r')
        ax5.axvline(x=minF4s,color= 'aqua')
        ax5.set_ylabel('$F4prim$')
        ax0.set_title('tr'+component+' Windows time '+ str(T)+'s')
        plt.savefig('/home/claire/PHD/Working/WavePicking_%s.png'%component,dpi=300)
    return  Second-60*5+minF4p/100, Second-60*5+minF4s/100
    
def WavePicking2(tr,T, SecondP,SecondS,plot):
    """  give start of p and s respecting the theorical time delay between P and S calculated using a Chi 2001 taiwanese velocity model (surface velocities) and Iasp91  (deep velocities)
    -input:
        -tr : trace 
        -T : float, 5s time of the sliding windows see baillard but could be change trial shoaw goiod result with 5s. 
        -plot: bool, if true plot will be draw
        -SecondP, float theorical P arrival in second 
        -Seconds,float theorical S arrival in second 
    -output
        - minF4p: float, coordiantes of the p start in second
        - minF4s: float, coordiantes of the s start in second
    -exemple:
        minF4p, minF4s = WavePicking(tr,5,False)"""
    tr_copy = tr.copy()
    tr_copy2 = tr.copy()
    P_Sdelay = (SecondS-SecondP)
  
    #NESURTOUT PAS FILTERER
    #tr_copy.filter('bandpass', freqmin=0.2, freqmax= 15)
    if SecondP > 60: #second 
        tr_copy.trim(starttime =tr.stats.starttime + SecondP-60)
#    if 3600 > Second+60*5:
#        tr_copy.trim(endtime =tr.stats.starttime + Second+60*10)
    
    if (tr.stats.starttime +SecondP+P_Sdelay+120)<tr.stats.endtime:
        print 'end time', (tr.stats.starttime +SecondP+P_Sdelay+120)
        tr_copy.trim(endtime =tr.stats.starttime + SecondP+P_Sdelay+120)

    Cf = CFkurtosis(tr_copy,T) 

    f2 = F2(Cf)

    f3 = F3(f2)

    f3prim = MeanF3(f3,0.5) #sliding windows of 0.5 s show good results 
 
    f4 = F4(f3prim)
   
    # P and S coordinates 
    minF4p = np.argmin(f4[0: len(f4)-100])
    if P_Sdelay <>0 :
        if len(f4[minF4p+int((1-0.20)*P_Sdelay*100):minF4p+int((1+0.50)*P_Sdelay*100)])<>0:
        # min s wave picking between after p plus 50% of the therorical  ps delay
            min2F4s2 = np.argmin(f4[minF4p+int((1-0.20)*P_Sdelay*100):minF4p+int((1+0.50)*P_Sdelay*100)])
            minF4s2 = minF4p+(1-0.20)*P_Sdelay*100+min2F4s2
        else : 
            minF4s2 = minF4p+10
    else :
        min2F4s2 = np.argmin(f4[minF4p+int((1-0.20)*2*100):minF4p+int((10)*2*100)])
        minF4s2 = minF4p+(1-0.20)*2*100+min2F4s2
    #tr_copy2.trim(starttime=tr_copy2.stats.starttime + SecondP-60 + int(minF4p/100)+ P_Sdelay)
  
    #minimum of points between p and s start function of distance  (calculated with vs =2.6 km.s-1 and vp = 4.1km.s-1)
 
    if plot == True:
        component = tr.stats.channel
        f, (ax0,ax1, ax2, ax3, ax4,ax5) = plt.subplots(6, sharex=True)
        x = range(len(tr_copy.data))
        x2 = range(len(tr_copy2.data))
        ax0.plot(x,tr_copy.data)
        ax0.axvline(x=minF4p,color= 'r')
        ax0.axvline(x=minF4s2,color= 'aqua')
        ax0.set_ylabel('$velocity m.s^-1$')
        ax1.plot(x,Cf)
        ax1.axvline(x=minF4p,color= 'r')
        ax1.axvline(x=minF4s2,color= 'aqua')
        ax1.set_ylabel('$Cf$')
        ax2.plot(x, f2)
        ax2.axvline(x=minF4p,color= 'r')
        ax2.axvline(x=minF4s2,color= 'aqua')
        ax2.set_ylabel('$F2$')
        ax3.plot(x,f3prim)
        ax3.set_ylabel('$F3$')
        ax3.axvline(x=minF4p,color= 'r')
        ax3.axvline(x=minF4s2,color= 'aqua')
        ax4.plot(x, f4)
        ax4.axvline(x=minF4p,color= 'r')
        ax4.axvline(x=minF4s2,color= 'aqua')
        ax4.set_ylabel('$F4$')
        ax0.set_title('tr'+component+' Windows time '+ str(T)+'s')
        plt.savefig('/home/claire/PHD/Working/WavePicking_%s.png'%component,dpi=300)
    return  SecondP-60+minF4p/100,SecondP-60+minF4s2/100