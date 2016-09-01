# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 23:43:48 2016

@author: claire
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


"""List of functions: Event selections
- SelectEventAz : select event from event caracterisation txt file according to their azimuths for a given station
- SelectEventML : select event from event caracterisation txt file according to their magnitude ML for a given station
- SelectEventDist : select event from event caracterisation txt file according to their radiale distance to stations for a given station
- SelectEventMLDist : Return an array containing the station -event caracteristics for events with a magnitude dbetween [MlMin; MlLMax] and distant from [DistMin, DistMax]
- SelectEventJJul : 
- ConvertDatestr: return year , station and hour in str format as use to select SAC original files 

"""

def SelectEventAZ(File,Azmin, Azmax, station):
    """Return an array containing inforlations about all events that are between Azmin and Azmax from the station
    * inputs : 
        - File: txt file containing event informations
        - Azmin type: float , min Azimuth
        - Azmax type  float , max Azimuth 
        - station: reference station
    * output :
       - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st of the event which have Azimuth betweeen [Azmin;Azmax]
        where: 
                - Station: station considered
                - year :  year
                - julian_day : Julian day 
                - Hour :  hour of the event
                - seconds : event start time of the P wave in second calculated with iasp model 
                - ML :  magnitude 
                - depth :  depth in km
                - Rev_st: radiale distance between the source and the event in km 
                - BAz_ev_st :  back azimuth
                - Az_ev_st :  azimuth 
        - Aevent.shape: type list , shape of the array 
        - column : column of Azimuth      
    * exemple: 
     File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
     SelectEventAZ(File,-179, -170, '7')
        ==> [[  7.00000000e+00   2.01500000e+03   1.00000000e+00 ...,   6.20355793e+01
                9.25351662e+00  -1.70710661e+02]
             [  7.00000000e+00   2.01500000e+03   2.00000000e+00 ...,   2.73034494e+02
                4.72828341e+00  -1.75191203e+02]
             [  7.00000000e+00   2.01500000e+03   2.00000000e+00 ...,   1.63196327e+02
                3.25224440e+00  -1.76712617e+02]
             ..., 
             [  7.00000000e+00   2.01500000e+03   3.65000000e+02 ...,   1.91452315e+02
                2.46917922e+00  -1.77499746e+02]
             [  7.00000000e+00   2.01500000e+03   3.65000000e+02 ...,   2.66364290e+02
                3.03479811e+00  -1.76915277e+02]
             [  7.00000000e+00   2.01500000e+03   3.65000000e+02 ...,   2.74222002e+02
                8.66057049e+00  -1.71186326e+02]], ((1425, 10))"""
    
    
     #1.  read file and find the column of intereste number~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B = open(File,'r')
    for ligne in B: 
        L1 = ligne.split(',')
        print L1
        columnAz = L1.index('Az\n')
        break
    B.close()
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
   
    #2. select elevent of the array between Mwmin and Mw max : CoordSelected2 [[lign, column]]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CoordSelected =  np.argwhere((A[:,columnAz]>=Azmin)&( A[:,columnAz]<Azmax)&( A[:,0]==float(station)))
     
    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent=[]
    Levent = [A[i] for i in CoordSelected[:,0] if A[i] not in Levent]
     
    Aevent = np.array(Levent)
    return Aevent, Aevent.shape, columnAz
       
def SelectEventML(File,Mwmin, Mwmax, station):
    """ Return an array containing inforlations about all events that havce a magnitude between MLmin and MLmax and its shape
    * inputs : 
        - File: txt file containing event informations
        - Azmin type: float , min Magnitude
        - Azmax type  float , max Magnitude
        - station: reference station
    * output :
        - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st of the event that are  betweeen [Mlmin , MLmax]
            where: 
                    - Station: station considered
                    - year :  year
                    - julian_day : Julian day 
                    - Hour :  hour of the event
                    - seconds : event start time of the P wave in second calculated with iasp model 
                    - ML :  magnitude 
                    - depth :  depth in km
                    - Rev_st: radiale distance between the source and the event in km 
                    - BAz_ev_st :  back azimuth
                    - Az_ev_st :  azimuth 
        - Aevent.shape: type list , shape of the array 
        - Cml = column of Ml 
    * exemple: 
    File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'    
    SelectEventML(File,3, 4, '7')"""
   
    #1.  read file and find the column of intereste number~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B = open(File,'r')

    for ligne in B: 
        L1 = ligne.split(',')
        columnMw = L1.index('ML')
        break
    B.close()
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
     
    #2. select elevent of the array between Mwmin and Mw max : CoordSelected2 [[lign, column]]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CoordSelected =  np.argwhere((A[:,columnMw]>=Mwmin)&( A[:,columnMw]<Mwmax)&( A[:,0]==float(station)))
    print "ll"
    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent = []
    Levent = [A[i] for i in CoordSelected[:,0] if A[i] not in Levent]
    
    Aevent = np.array(Levent)
    return Aevent, Aevent.shape, columnMw

def SelectEventDist(File,Distmin, Distmax, station):
    """ Return an array containing inforlations about all events that are distante from Distmin to Distmax of the station 
        - File: txt file containing event informations
        - Distmin type: float , min radiale Distance between station and event 
        - Distmax type  float , max radiale Distance between station and event 
        - station: reference station
    * output :
        - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st lat lobng of the event that are  betweeen [Distmin ,Distmax]
            where: 
                    - Station: station considered
                    - year :  year
                    - julian_day : Julian day 
                    - Hour :  hour of the event
                    - seconds : event start time of the P wave in second calculated with iasp model 
                    - ML :  magnitude 
                    - depth :  depth in km
                    - Rev_st: radiale distance between the source and the event in km 
                    - BAz_ev_st :  back azimuth
                    - Az_ev_st :  azimuth 
		    - Lat; latitude of event
		    - Long : longitude
        - Aevent.shape: type list , shape of the array 
        -  columnDist :  type int : column of the distance
    * exemple: 
        File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
        SelectEventDist(File,10, 20, '7')"""
   
    #1.  read file and find the column of intereste number~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B = open(File,'r')
    for ligne in B: 
       L1 = ligne.split(',')
       columnDist = L1.index('R')
       break
    B.close()
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
       
    #2. select elevent of the array between Mwmin and Mw max : CoordSelected2 [[lign, column]]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CoordSelected =  np.argwhere((A[:,columnDist]>=Distmin)&( A[:,columnDist]<Distmax)&( A[:,0]==float(station)))
     
    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent=[]
    Levent = [A[i] for i in CoordSelected[:,0] if A[i] not in Levent]
     
    Aevent = np.array(Levent)
    return Aevent, Aevent.shape, columnDist
    

def SelectEvent(File,station):
    """ Return an array containing the station 
    * input : 
        - File: txt file containing event informations
        - station: reference station
    * output :
        - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st where we have satas about a given station
            where: 
                    - Station: station considered
                    - year :  year
                    - julian_day : Julian day 
                    - Hour :  hour of the event
                    - seconds : event start time of the P wave in second calculated with iasp model 
                    - ML :  magnitude 
                    - depth :  depth in km
                    - Rev_st: radiale distance between the source and the event in km 
                    - BAz_ev_st :  back azimuth
                    - Az_ev_st :  azimuth 
		    - Lat; latitude of event
		    -Long : longitude
        - Aevent.shape: type list , shape of the array 
       
    * exemple: 
        File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
        SelectEventDist(File,10, 20, '7')"""
      
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
       
    #2. select ligne corresponding to the given station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CoordSelected =  np.argwhere(A[:,0]==float(station))
     
    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent=[]
    Levent = [A[i] for i in CoordSelected[:,0] if A[i] not in Levent]
     
    Aevent = np.array(Levent)
    return Aevent, Aevent.shape

def SelectEventMLDist(File, station, MlMin, MlMax, DistMin, DistMax) : 
    """ Return an array containing the station -event caracteristics for events with a magnitude dbetween [MlMin; MlLMax] and distant from [DistMin, DistMax]
     
    * input : 
        - File: txt file containing event informations
        - station: reference station
        -MlMin type float, Minimum magnitude
        - MlMax, type Float, Maximum Magnitude
        -DistMin, Type Minimum distance 
        -DistMax , type float Maximumdistance 
        -output, type str; file where the final aray ll be save 
    * output :
        - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st where we have satas about a given station
            where: 
                    - Station: station considered
                    - year :  year
                    - julian_day : Julian day 
                    - Hour :  hour of the event
                    - seconds : event start time of the P wave in second calculated with iasp model 
                    - ML :  magnitude 
                    - depth :  depth in km
                    - Rev_st: radiale distance between the source and the event in km 
                    - BAz_ev_st :  back azimuth
                    - Az_ev_st :  azimuth 
		    - Lat; latitude of event
		    -Long : longitude
        - Aevent.shape: type list , shape of the array 
       
    * exemple: 
        File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
        SelectEventDist(File,10, 20, '7')"""
        
    List=[]
    B = open(File,'r')
    for ligne in B: 
       L1 = ligne.split(',')
       columnDist = L1.index('R')
       columnMl = L1.index('ML')
       break
    for i, line in enumerate(B):
        L1=line.split(',')
        if 'a' in L1[0]:
            List.append(i)
        
    B.close()
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
       
    #2. select ligne corresponding to the given station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    CoordSelectedinf4 =  np.argwhere((A[:,0]==float(station)) & (A[:,0]<>'5a') & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<100) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<4))
#    CoordSelectedsup4 =  np.argwhere((A[:,0]==float(station)) & (A[:,0]<>'5a') & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=4) & (A[:, columnMl]<MlMax))
    CoordSelected =  np.argwhere((A[:,0]==float(station)) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<MlMax))
    #CoordSelected =np.concatenate((CoordSelectedinf4,CoordSelectedsup4)) 

    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent=[]
    Levent = [A[i] for i in CoordSelected[:,0] if A[i] not in Levent]
     
    Aevent = np.array(Levent)
   
    return Aevent, Aevent.shape
    
def SelectEventMLDistSt3(File, station, MlMin, MlMax, DistMin, DistMax) : 
    """ Return an array containing the station -event caracteristics for events with a magnitude dbetween [MlMin; MlLMax] and distant from [DistMin, DistMax] to station 3
     
    * input : 
        - File: txt file containing event informations
        - station: reference station
        -MlMin type float, Minimum magnitude
        - MlMax, type Float, Maximum Magnitude
        -DistMin, Type Minimum distance 
        -DistMax , type float Maximumdistance 
        -output, type str; file where the final aray ll be save 
    * output :
        - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st where we have satas about a given station
            where: 
                    - Station: station considered
                    - year :  year
                    - julian_day : Julian day 
                    - Hour :  hour of the event
                    - seconds : event start time of the P wave in second calculated with iasp model 
                    - ML :  magnitude 
                    - depth :  depth in km
                    - Rev_st: radiale distance between the source and the event in km 
                    - BAz_ev_st :  back azimuth
                    - Az_ev_st :  azimuth 
		    - Lat; latitude of event
		    -Long : longitude
        - Aevent.shape: type list , shape of the array 
       
    * exemple: 
        File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
        SelectEventDist(File,10, 20, '7')"""
    Dictst = {'1':-2,'2':-1,'3':0,'4':1,'5':2,'6':3,'7':4}   
    List=[]
    B = open(File,'r')
    for ligne in B: 
       L1 = ligne.split(',')
       columnDist = L1.index('R')
       columnMl = L1.index('ML')
       break
    for i, line in enumerate(B):
        L1=line.split(',')
        if 'a' in L1[0]:
            List.append(i)
        
    B.close()
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
       
    #2. select ligne corresponding to the given station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    CoordSelectedinf4 =  np.argwhere((A[:,0]==float(station)) & (A[:,0]<>'5a') & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<100) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<4))
#    CoordSelectedsup4 =  np.argwhere((A[:,0]==float(station)) & (A[:,0]<>'5a') & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=4) & (A[:, columnMl]<MlMax))
    CoordSelected =  np.argwhere((A[:,0]==3) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<MlMax))
    CoordSelectedst = CoordSelected + Dictst[station]
    print 'coord3',len(CoordSelected), 'coordst',len(CoordSelectedst)
    #CoordSelected =np.concatenate((CoordSelectedinf4,CoordSelectedsup4)) 

    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent=[]
    Levent = [A[i] for i in CoordSelectedst[:,0] if A[i] not in Levent]
     
    Aevent = np.array(Levent)
     
    return Aevent, Aevent.shape
    
    
    
def Array_St_EQ_MlDist(File, MlMin, MlMax, DistMin, DistMax) : 
    """ Return an array containing the station -event caracteristics for events with a magnitude dbetween [MlMin; MlLMax] and distant from [DistMin, DistMax] to station 3
     
    * input : 
        - File: txt file containing event informations
        - station: reference station
        -MlMin type float, Minimum magnitude
        - MlMax, type Float, Maximum Magnitude
        -DistMin, Type Minimum distance 
        -DistMax , type float Maximumdistance 
        -output, type str; file where the final aray ll be save 
    * output :
        - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st where we have satas about a given station
            where: 
                    - Station: station considered
                    - year :  year
                    - julian_day : Julian day 
                    - Hour :  hour of the event
                    - seconds : event start time of the P wave in second calculated with iasp model 
                    - ML :  magnitude 
                    - depth :  depth in km
                    - Rev_st: radiale distance between the source and the event in km 
                    - BAz_ev_st :  back azimuth
                    - Az_ev_st :  azimuth 
		    - Lat; latitude of event
		    -Long : longitude
        - Aevent.shape: type list , shape of the array 
       
    * exemple: 
        File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
        SelectEventDist(File,10, 20, '7')"""
    Dictst = {'1':-2,'2':-1,'3':0,'4':1,'5':2,'6':3,'7':4}   
    List=[]
    B = open(File,'r')
    for ligne in B: 
       L1 = ligne.split(',')
       columnDist = L1.index('R')
       columnMl = L1.index('ML')
       break
    for i, line in enumerate(B):
        L1=line.split(',')
        if 'a' in L1[0]:
            List.append(i)
        
    B.close()
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})      
    #2. select ligne corresponding to the given station first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CoordSelected =  np.argwhere((A[:,0]==3) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<MlMax))
    #take the coordinates of the EQ corresponding to the station considered
    CoordSelectedst = CoordSelected
    for station in ['1','2','4','5','6','7']:
        CoordSelectedst = np.concatenate((CoordSelectedst,CoordSelected + Dictst[station]))
        
    print 'coord3',len(CoordSelected), 'coordst',len(CoordSelectedst)
    #CoordSelected =np.concatenate((CoordSelectedinf4,CoordSelectedsup4)) 

    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent=[]
    Levent = [A[i] for i in CoordSelectedst[:,0] if A[i] not in Levent]
     
    Aevent = np.array(Levent)
    LEQ= [str()]
   
    return Aevent

def DeleteEventMLDistJJulref(File, station, MlMin, MlMax, DistMin, DistMax):
    """ Return an array containing the station -event caracteristics for events with a magnitude dbetween [MlMin; MlLMax] and distant from [DistMin, DistMax] for a given julian day reference 
     
    * input : 
        - File: txt file containing event informations
        - station: reference station
        -MlMin type float, Minimum magnitude
        - MlMax, type Float, Maximum Magnitude
        -DistMin, Type Minimum distance 
        -DistMax , type float Maximumdistance 
        -output, type str; file where the final aray ll be save 
    * output :
        - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st where we have satas about a given station
            where: 
                    - Station: station considered
                    - year :  year
                    - julian_day : Julian day 
                    - Hour :  hour of the event
                    - seconds : event start time of the P wave in second calculated with iasp model 
                    - ML :  magnitude 
                    - depth :  depth in km
                    - Rev_st: radiale distance between the source and the event in km 
                    - BAz_ev_st :  back azimuth
                    - Az_ev_st :  azimuth 
		    - Lat; latitude of event
		    -Long : longitude
        - Aevent.shape: type list , shape of the array 
       
    * exemple: 
        File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
        SelectEventDist(File,10, 20, '7')"""
    List=[]
    B = open(File,'r')
    for ligne in B: 
       L1 = ligne.split(',')
       columnDist = L1.index('R')
       columnMl = L1.index('ML')
       break
    for i, line in enumerate(B):
        L1=line.split(',')
        if 'a' in L1[0]:
            List.append(i)
        
    B.close()
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
    A =np.delete(A, List, axis=0)
    #2. select ligne corresponding to the given station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # CoordSelectedinf4 =  np.argwhere((A[:,0]==float(station))  & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<50) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<4) )
    #CoordSelectedMid4 =  np.argwhere((A[:,0]==float(station)) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<100) & (A[:, columnMl]>=4) & (A[:, columnMl]<5) )
    #CoordSelectedsup4 =  np.argwhere((A[:,0]==float(station))  & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=5) & (A[:, columnMl]<MlMax))
    #2. select ligne corresponding to the given station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    CoordSelectedinf4 =  np.argwhere((A[:,0]==float(station)) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<50) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<4))
#    CoordSelectedMid4 =  np.argwhere((A[:,0]==float(station)) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<100) & (A[:, columnMl]>=4) & (A[:, columnMl]<5))
#    CoordSelectedsup4 =  np.argwhere((A[:,0]==float(station)) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=5) & (A[:, columnMl]<MlMax))
    #CoordSelectedsup4 =  np.argwhere((A[:,0]==float(station)) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=5) & (A[:, columnMl]<MlMax))
    #CoordSelected =np.concatenate((CoordSelectedinf4,CoordSelectedsup4,CoordSelectedMid4)) 
    CoordSelected =  np.argwhere((A[:,0]==float(station)) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<MlMax))
    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent=[]
    Levent = [A[i] for i in CoordSelected[:,0] if A[i] not in Levent]
    CoordSelected =np.concatenate((CoordSelectedinf4,CoordSelectedsup4,CoordSelectedMid4)) 
    LJul = [A[i,2] for i in CoordSelected]
    for i in xrange(len(A)):
        if A[i,2] in LJul : 
            A = np.delete(A,i,axis=0) #delete all ligne
        
    
    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    return A, A.shape

def SelectEventMLDistJJulref(File, station, MlMin, MlMax, DistMin, DistMax, LJulref):
    """ Return an array containing the station -event caracteristics for events with a magnitude dbetween [MlMin; MlLMax] and distant from [DistMin, DistMax] for a given julian day reference 
     
    * input : 
        - File: txt file containing event informations
        - station: reference station
        -MlMin type float, Minimum magnitude
        - MlMax, type Float, Maximum Magnitude
        -DistMin, Type Minimum distance 
        -DistMax , type float Maximumdistance 
        -output, type str; file where the final aray ll be save 
    * output :
        - Aevent :  type np array each line gives  : station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st lat, long where we have satas about a given station
            where: 
                    - Station: station considered
                    - year :  year
                    - julian_day : Julian day 
                    - Hour :  hour of the event
                    - seconds : event start time of the P wave in second calculated with iasp model 
                    - ML :  magnitude 
                    - depth :  depth in km
                    - Rev_st: radiale distance between the source and the event in km 
                    - BAz_ev_st :  back azimuth
                    - Az_ev_st :  azimuth 
        - Aevent.shape: type list , shape of the array 
       
    * exemple: 
        File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
        SelectEventDist(File,10, 20, '7')"""
    List=[]
    B = open(File,'r')
    for ligne in B: 
       L1 = ligne.split(',')
       columnDist = L1.index('R')
       columnMl = L1.index('ML')
       break
    for i, line in enumerate(B):
        L1=line.split(',')
        if 'a' in L1[0]:
            List.append(i)
        
    B.close()
    A = np.loadtxt(File,skiprows=1,delimiter = ',',converters ={0:lambda x: int(x[1])})
    A =np.delete(A, List, axis=0)
    #2. select ligne corresponding to the given station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print type(A[:,2]),A[:,2],A[:,1]
    CoordSelectedinf4 =  np.argwhere((A[:,0]==float(station))  & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<50) & (A[:, columnMl]>=MlMin) & (A[:, columnMl]<4) & (A[:,2] == float(LJulref)))
    CoordSelectedMid4 =  np.argwhere((A[:,0]==float(station)) & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<100) & (A[:, columnMl]>=4) & (A[:, columnMl]<5) & (A[:,2] == float(LJulref)))
    CoordSelectedsup4 =  np.argwhere((A[:,0]==float(station))  & (A[:, columnDist]>=DistMin) & (A[:, columnDist]<DistMax) & (A[:, columnMl]>=5) & (A[:, columnMl]<MlMax) & (A[:,2] ==float(LJulref)))
   
    CoordSelected =np.concatenate((CoordSelectedinf4,CoordSelectedsup4,CoordSelectedMid4)) 

    #3. list of events :  build a new array containing all the information about the event we are interrested in~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Levent=[]
    Levent = [A[i] for i in CoordSelected[:,0] if A[i] not in Levent]
     
    Aevent = np.array(Levent)
   
    return Aevent, Aevent.shape
    

def ConvertDatestr(ligne):
    """ return year , station and hour in str format as use to select SAC original files 
    * input:  
        - ligne :  ligne of an array after selection or without selection ie 1D array
    * output :
        - year : type str; year of the event consider ex : 2015e+3 ==>  15 
        - station : type str ; station considered ex 1e+3 ==> 1
        - Hour type str; hour of the event ex 1 ==> 01   
    * exemple :
        ConvertDatestr([1.00000000e+00, 2.01500000e+03, 1.00000000e+00, 0.00000000e+00,3.21009088e+02,1.60000000e+00,9.23000000e+00,6.97429071e+01,-1.75018665e+02,4.95750780e+00])
        ==>('1', '15', '01')
    """
    print ligne
    station = str(int(ligne[0]))
    year = str(int(ligne[1])).split('0')[1]
    
    Second =int(ligne[4])
    ml = float(ligne[5])
    depth = float(ligne[6])
    Rdistance = float(ligne[7])
    Lat = float(ligne[10])
    Long = float(ligne[11])
    Az = float(ligne[9])
    if Second>3600 : 
        h = int(Second)//3600
        Second = int(Second)%3600
        hour = int(ligne[3])+h 
        ml = float(ligne[5])
        if hour<10:
            hour ='0'+str(hour)
        else :
            hour = str(hour)
    else : 
        if ligne[3]<10:
            hour ='0'+str(int(ligne[3]))
        else :
            hour = str(int(ligne[3]))
    if ligne[2]<10:
        jJul = '00'+str(int(ligne[2]))
    elif ligne[2]<100:
        jJul = '0'+str(int(ligne[2]))
    else :
        jJul = str(int(ligne[2]))
    return station, year,jJul, hour, Second,ml,depth, Rdistance, Lat,Long, Az




#TEts ################################################################
#B= open(File,'r')
#i=0

#for ligne in B :
#    if i<= 0: 
#        L = ligne.split('/n')
#        print len(ligne),len(L)
#        LB = np.array(L)
#        print 'done', type(LB) ,LB[0], LB[1], LB[7]
#        np.savetxt('/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt',LB,fmt = '%s')
#        print 'lll'
#        i+=1
#    else :
#        break
#B.close()
#output = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics_ML%s_%s_Dist_%s_%s.txt'%(3, 9, 0, 200)
#File ='/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
#print SelectEventMLDist(File, '1',3, 9, 0, 200, output)
#print ConvertDatestr([1.00000000e+00, 2.01500000e+03, 1.00000000e+00, 0.00000000e+00,3.21009088e+02,1.60000000e+00,9.23000000e+00,6.97429071e+01,-1.75018665e+02,4.95750780e+00])
