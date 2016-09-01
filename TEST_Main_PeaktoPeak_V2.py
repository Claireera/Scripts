# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 08:09:27 2016

@author: claire

New Peak to peak test: This code is use to show the diffrent step to obtain the peak to peak of a given signal 
  
"""

#from Events_caracterisation import *
from Events_Selections import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *
from Instance_EQPeaktoPeak import *
from Signal_Wave_Picking import *
import cPickle
import json


#Earthquake features 
MlMin= 4
MlMax = 5
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'
LStations = ['1','2','3','4','5','6','7']
#component 
LComponents = ['N','R','E','T','Z','H']
# 0. limit SNR above signal is valid
LimSNR = 2
with open('/home/claire/PHD/Working/Data/List_EQ_St_Comp_Freq_Ml_%s_%s_distmax_%s_Test.txt'%(MlMin,MlMax,DistMax)) as f:
    LEqStCompFreq =json.load(f)
    
for EQstCompFreq in LEqStCompFreq[0:1][:]:
        #            # 6.0 take case feature EQ, station component and frequency
        EQname, Station, Year,jJul, Hour, Secondp, Seconds,Ml,Depth, Rdistance, Lat,Long, Az,BAz, Component, number, freqband = EQstCompFreq
             
            
#        if Year<>f[EQname].attrs['Year'] or  Hour<>f[EQname].attrs['Hour'] or  jJul<>f[EQname].attrs['JJul']: 
#            print 'ERROR PROBLEM OF EQ FEATURES' , 'Je suis le worker {} et je traitle le seisme {}, et la composante {} pour la frequence max {}'.format(me,EQname,Component, freqband[1])
            
            
        ## 6.1 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if not  os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
           print 'file not existing Year :',  Year, 'Julian day : ',jJul, 'Hour : ', Hour
           continue
       
        st=Read_event(Year,jJul,Hour,Secondp, Station, False)
        print Hour, type(Hour)
        #if the begining of the signal is too close to 00min then the wave peaking won't be accurate thus it's necessary to merge the signal whith the previous one! 
        if Secondp<100:
            if int(Hour)-1>= 0:
                st = st + Read_event(Year,jJul,str(int(Hour)-1),Secondp, Station, False)
            else : 
                if int(jJul)-1>1:
                    st = st + Read_event(Year,int(jJul-1),str(23),Secondp, Station, False)
                else : 
                    st = st + Read_event(Year-1,int(365),str(23),Secondp, Station, False)
        st.merge(method= 1)
        
        #6.2 Instrument correction~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        stcorrect = Stream_Correction(st, Station, False)
        
        #6.3 Rotation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        strot = Rotation(stcorrect,BAz, False)
        st = stcorrect.copy()

#           print 'There is still 
        #6.4 select the trace of the component~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if Component in ['R', 'T']:
            st_copy = st.copy()
            strot_copy= strot.copy()
            tr = strot_copy.select(component=Component)[0]
            trZ = st_copy.select(component="Z")[0]
            trE = st_copy.select(component="E")[0]
        elif Component=='H':
            st_copy = st.copy()
            strot_copy= strot.copy()
            tr, trZabs = RealafterrottrH(strot_copy,False)
            trZ = st_copy.select(component="Z")[0]
            trE = st_copy.select(component="E")[0]
        else :
            st_copy = st.copy()
            tr = st_copy.select(component=Component)[0]
            trZ = st_copy.select(component="Z")[0]
            trE = st_copy.select(component="E")[0]

        #6.5Define Python instance nameEQ,magnitude,depth,Rdistance,Lat,Long,Az,time,Station, frequency band, stream, streamrot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        EQ = PtPearthquake(EQname,Ml,Depth,Rdistance,Lat,Long,Az,[Year,jJul,Hour,str(int(Secondp)),str(int(Seconds))],Station,freqband,tr,Component)         
      
        #6.6 P and S arrivals calculated on the Z componenent~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print int(EQ.time[3]),int(EQ.time[4])
        EQ.second_P, second_S = WavePicking2(trZ,5,int(EQ.time[3]),int(EQ.time[4]),True)
        second_P, EQ.second_S = WavePicking2(trE,5,int(EQ.time[3]),int(EQ.time[4]),True)
        print int(EQ.time[3]),int(EQ.time[4]),EQ.second_S,EQ.second_P
        #6.7 Trim the signal between ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        EQ.trTrim =  Trimtr(EQ.tr,EQ.second_P,1.5, 0.65,7,15,False,False)

        #6.8SNR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        EQ.valid,EQ.SNR= SNRstd(EQ.tr, EQ.trTrim, LimSNR, EQ.second_P, False)
        #6.9saturation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #6.9bisplot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fig = plt.figure()
        if EQ.valid ==True:
            val = "T"
        else : 
            val ="F"
        
        #6.10 PGA, PGV, Arias, Peak to Peak~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if EQ.valid==True:
            #PGAG
            EQ.FPGA()
            #PGV
            EQ.FPGV()
            #Arias
            EQ.FArias()
            #2.6 .4 Peak to Peak and filter
            EQ.PtP_envelope_Test()
            
        
        plt.figure()
        plt.plot(EQ.trTrim.data, color='k')
        plt.plot(EQ.trFilt.data,color = 'b',linestyle ='--')
        plt.axvline(x=0,color= 'orchid')
        plt.axvline(x=(EQ.second_S-EQ.second_P)*100,color= 'aqua')
        plt.title (' valid'+val+'Start '+str(EQ.second_P)+ 'EQ '+ EQname  +'component ' + Component + 'st '+Station + ' peaktopeak:'+ str(EQ.PeaktoPeak))
        plt.savefig('/home/claire/PHD/Working/Draft/EQ_'+ EQ.nameEQ  +'_component_' + EQ.Component + '_st_'+EQ.Station+'Fmax_'+str(EQ.frband[1])+'ML_'+str(MlMax)+'.png')
        print 'Snr calcul',EQ.SNR, val, 'valid'
         
#                2.7 Dictionary containing parameters definition of the EQ considered at the given station for a given frequency of filtering
#                dataPickle = qstatcPickle.dumps(EQ)
#        Dict = EQ.__dict__
#
f.close()
        
