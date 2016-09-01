# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 15:42:13 2016

@author: claire
"""


import numpy as np
from math import radians
from pyproj import Geod
from obspy.taup import TauPyModel
from obspy.taup.taup import getTravelTimes
from obspy.core.utcdatetime import UTCDateTime

"""" list of functions : Event caracterisations
- Epidist : returns surface distance in km, and in degree and the real distance between the source and the station considered 
- Azimuth : returns the azimuth back azimuth and distance from the source to the event
- jJulian : Jour julian, hour and second of a date yyyymmddhhmmss
- ParrivalTime(Depthev,distangular) : return P arrival time
- Catalogue_file(FileStation, inputFile, outputfile) :  write an output file containing informations between station and events

"""

pi = 3.141592653589793
def Epidist(latsource,longsource,latstation,longstation,depth):
    """distance between the point 1 and point 2 in a cartesian system 
    * inputs :
        -latsource, longsource : type float;  latitude and longitude of the event
        -latstation, longsation:  type float;  latitude and longitude of the station
        -depth:  type float; depth of the event in km 
    * output : 
        - Distkm :  type float; Surfacial distance between the source and the event in km
        - Distangle : type float; Surfacial distance between the source and the event in degree
        - Realdistance :  type float; radial distance between the station and the event (depth is taken into account)
        
    * exemple : 
        Epidist(22.3, 121.1, 23, 122, 12)
        ==>(126.78172275582767, 1.140175425099142, 127.34836168924808)
        """
    Distangle = np.sqrt((latsource-latstation)**2+(longsource-longstation)**2) 
    distangle = np.arccos(np.sin(latsource)*np.sin(latstation)+np.cos(latsource)*np.cos(latstation)*np.cos(latsource-latstation))
    Realdistance =  np.sqrt(((latsource-latstation)*(pi/180)*6371)**2+((longsource-longstation)*(pi/180)*6371)**2+depth**2) 
    DistKm = (pi/180)*Distangle*6371
    return DistKm, Distangle, Realdistance, distangle



def Azimuth(latsource, longsource, latsat,longsat, depth):
    """returns backAzimuth , azimuth in degrees [-180,180], and distance in km beween event and station
    * inputs :
        - latsource, longsource :type : float ; latitude and longitude of the event
        - latstation, longsation: type : float ; latitude and longitude of the station
        - depth: type : float ; depth of the event in km 
    * output :
        - BA: type : float ; Aack azimuth between the source and the station
        - A :type : float ; Azimuth azimuth between the source and the station
        - Dist : type : float ; surfaciam distance between the source and the station in m
        - Realdistance : type : float ; radial distance between the station and the event (depth is taken into account)
        - dcartangletype : float ; real distance in 
    * exemple : 
        Azimuth(22.3, 121.1, 23, 122, 12)
        ==>(49.8655392179202, -129.78785743989272, 120690.53790617369, 127.34836168924808)
        """
   
    dcartkm, dcartangle,Realdistance,distangle = Epidist(latsource,longsource,latsat,longsat,depth)
    g = Geod(ellps='WGS84')
    #surfacial distance, azimuth and backazimuth in degree between station and source
    BA,A,Dist = g.inv(longsource,latsource,longsat,latsat,radians=False) # BA back azimut (Event,station),A: azimut (station,Event) 
    return BA, A, Dist, Realdistance,dcartangle,distangle
    
    
def jJulian(date):
    """ returns year, jjul, hour in UTC, event Start time in second
    * inputs : 
        - date: written as yyyymmddhhmmss
    * outputs : 
        - year : type : int ; year 
        - jjul : type : int ; hour
        - hour : type : int hour 
        - Start : type int; Second + minute *60
    * exemple :
        jJulian('20151111124052')
        ==> (2015, 315, 4, 2452)"""
    
    dt= UTCDateTime(date)-(3600*8)
    jjul =  dt.julday
    hour =  dt.hour
    #strat event in second in UTC 
    Start = dt.minute *60 + dt.second
    year=dt.year
    return year, jjul, hour, Start
  


  
def ParrivalTime(Depthev,distangular):
    """calculate the P arrival time with the iasp91 model
    * inputs : 
        - Depthev :  depth of event in km
        - distangular :  distance angular between source and event in degrees
    * outputs : 
        - a :  P arrival time (P wave time travel between source and earthquake)
    * exemple :
    
        ParrivalTime(10,80)
        ==> 729.55308823 % 729 second from the source to arrive at the station
    """
    
    model = TauPyModel(model="taiwan_chietal2001")                          
    arrivals= model.get_travel_times(source_depth_in_km= Depthev , distance_in_degree=distangular,phase_list = ["p","s"])
    if len(arrivals)==0 :
        P=0
        S=0
    else :
        if len(arrivals)==2: 
            P=arrivals[0].time
            S= arrivals[1].time
        else :
            #if s not detected S will arrive a 25ms after the P wave
            P=arrivals[0].time
            S= arrivals[0].time + 0.5
    return P,S #arrival of the P wave in second



def Catalogue_file(FileStation, inputFile, outputfile):
    """write a txt file containing 
        - Station: station considered
        - year :  year
        - julian_day : Julian day 
        - Hour :  hour of the event
        - seconds : event start time of the P wave in second calculated with iasp model 
        - ML :  magnitude 
        - depth :  depth in km
        - R : radiale distance between the source and the event in km 
        - BAZ :  back azimuth
        - Az :  azimuth 
    * inputs : 
        - FileSation :  txt file containing informations relative to the station (name, latitude and longitude)
        - inputfile : txt file containing informations relative to the vent (date, lat, long, depth, Ml)
        - Outputfile :txt file where the information are goning to be recorded 
        
    * outputs: 
        - File containing per ligne station year jJul hour Start ML depth Rev_st BAz_ev_st Az_ev_st
    * exemple : 
        FileStation = '/home/claire/PHD/Working/Data/station_infos.txt'
        inputFile ='/home/claire/PHD/Working/Data/2015Taiwan_catalogue.txt'
        outputfile = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
        Catalogue_file(FileStation, inputFile, outputfile)"""
        
    Inputfile =open(inputFile, 'r')
    Catstation=open(FileStation, 'r')
    outfile= open(outputfile,'wb')
    outfile.write("Station,year,julian_day, Hour,Start,ML,depth,R,BAZ,Az, LAte, Longe\n")
    #read file ligne by ligne   
    k=0
    for ligne in Inputfile : 
        Listcaract = ligne.split('\t')
        if k>0 : 
            #event informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            year, jJul, hour, start = jJulian(Listcaract[0])
          
            longevnt=float(Listcaract[1])
            latevent= float(Listcaract[2])
           
            #event depth in km
            depth = float(Listcaract[3])
            #magnitude        
            ML = float(Listcaract[4].split('\r')[0])
            b= 0 
            
            Catstation=open(FileStation, 'r')
            
            for LigneSt in Catstation:
                    #print LigneSt
                    FLigneSt = LigneSt.split(' ')
                   
                    #Station informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 
                    station = FLigneSt[0].split('W')[1]
                    
                    Latst = float(FLigneSt[4])
                   
                    Longst = float(FLigneSt[8])
                    
                    #Backazimuth, Azzimuth and distance between event and station, radiale distance between sation and event~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    BAev_st, Az_ev_st, Dist, Rev_st,Distangle,distangle =  Azimuth(latevent, longevnt, Latst, Longst, depth)
                    #Arrival of P waves ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    ParrTime,SarrTime = ParrivalTime(depth,distangle)
                   
                    Startp = start + ParrTime 
                    Starts = start + SarrTime 
                  
                    #Write in output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    outfile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(str(station),str(year),str(jJul), str(hour),str(Startp),str(Starts),str(ML),str(depth),str(Rev_st),str(BAev_st),str(Az_ev_st),str(latevent),str(longevnt)))
                 
            Catstation.close()
        k+=1
    
    outfile.close()
    Inputfile.close()
  
    return
#
#print Epidist(22.3, 121.1, 23, 122, 12)
#FileStation = '/home/claire/PHD/Working/Data/station_infos.txt'
#inputFile ='/home/claire/PHD/Working/Data/EQcatalog_2015_2016full.txt'
#outputfile = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics_2015_2016Full2.txt'
#Catalogue_file(FileStation, inputFile, outputfile)
#        