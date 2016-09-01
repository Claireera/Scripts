# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 21:52:40 2016

@author: claire
"""
""" Main Peak to Peak on moving window 
1. peak to peak ware calculate a save for a selection of events and then save in a file for each station (one ligne per event each column correcpond toa given band pass)
2. Ratio are made compare to a reference station typicaly '5' or '7'
3. Mean of the ratio are calculated and plot
 """

from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *

Dcolor={'1':'#2ba70f','2':'#3342ff','3':'#ff3333','4':'#ffac33','5':'#D7DF01','6':'#4B0082','7':'#FF00FF'}

MlMin= 3 
MlMax = 9
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'
LStations = ['1','2','3','4','5','6','7']
##########################################################################
#Phase 1 : Peak to Peak and envelope calculation
##########################################################################

File = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics_TEST.txt'
output = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics_Dist_%s_%s_Ml_%s_%s.txt'
for Station in LStations : 
    outfilePtpH = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_movingWindows_2Hz_0.2p_%s_H.txt'%Station
    outfilePEnvH = '/home/claire/PHD/Working/Results/Envelope/Envelopes_files/EnvelopeMax_movingWindows_2Hz_0.2p_%s_H.txt'%Station
    
    outfilePtpZ = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_movingWindows_2Hz_0.2p_%s_H_Z.txt'%Station
    outfilePEnvZ = '/home/claire/PHD/Working/Results/Envelope/Envelopes_files/EnvelopeMax_movingWindows_2Hz_0.2p_%s_H_Z.txt'%Station
    ### 1.0 select events following the chosen criterium (Mlmin and Mlmax and distance to the station)~~~~~~~~~~~~~~~~~~~~
    EventCaractSelec, shape = SelectEventMLDist(File, Station, MlMin, MlMax, DistMin, DistMax, output)
    print 'length event caract' , len(EventCaractSelec)
   # initialisation of list of list of Peak to Peak and event each line correspond to one event each column to a value obtain for a specific pass band filter
    AeventPtP=[]
    AeventEvelope=[]
    
    ###2. take event per event 
    for j in xrange(len(EventCaractSelec)):
     
        
        ###2.1 Event characteristics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Station, Year,jJul, Hour , Second,ml = ConvertDatestr(EventCaractSelec[j])
       BAz = EventCaractSelec[j][8]
       print  jJul, Hour , Year, 'year'
       
       ### 2.2 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
       if not  os.path.exists('/home/claire/PHD/Working/Data/Wanrong_array_2015/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
           print 'file not existing Year :',  Year, 'Julian day : ',jJul, 'Hour : ', Hour
           continue
       
       # plot control 
       if jJul == jJultoplot: 
           plot =True
       else : 
           plot= False
           
       st=Read_event(Year,jJul,Hour,Second, Station, plot)
       #change date jjul in 2016
       if int(Year)==16:
          jJul = int(jJul)+365
       else : 
          jJul= int(jJul)
          
       ### 2.3 instrumental correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       stcorrect = Stream_Correction(st, Station, False)
       del st
       
       ### 2.4 Rotation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print 'rotation'       
       strot = Rotation(stcorrect,BAz, plot)
       
       ### 2.6 filter on moving windows ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Fmini =0.1
       Fmaxi=2.
       
       #initialisation of Lptp and L eveloppe~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       LptpH = [jJul,Hour,ml]
       LenvmaxH = [jJul,Hour,ml]
       
       LptpZ = [jJul,Hour,ml]
       LenvmaxZ = [jJul,Hour,ml]
       
       for i in xrange(91):
           st_copy= strot.copy()
           Fmin= Fmini +0.2*i
           Fmax= Fmaxi +0.2*i
           # filter traces with moving windows 
           stSumtracefilt=st_copy.filter('bandpass',freqmin=Fmin,freqmax=Fmax,corners=4,zerophase=True)
                      
           #extract H quadratic average and  Z absolute components 
           trH, trZabs = RealafterrottrH(stSumtracefilt,plot)
           ### 2.7 Peak to Peak calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           #Horizontal 
           ptpmaxH,ptpFirstH= PeaktoPeak_positive(trH, plot)
           LptpH.append(ptpmaxH)
           #Vertical 
           ptpmaxZ,ptpFirstZ= PeaktoPeak_positive(trZabs, plot)
           LptpZ.append(ptpmaxZ)
           
           ### 2.8 Enveloppe ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           #Horizontal
           EnvH, tH, SumEnvelope, EnvMaxH = Envelope(trH,plot)
           LenvmaxH.append(EnvMaxH)
           #Vertical 
           EnvZ, tZ, SumEnvelopeZ, EnvMaxZ = Envelope(trZabs,plot)
           LenvmaxZ.append(EnvMaxZ)
       #horizontal
       AeventPtPH.append(LptpH)
       AeventEvelopeH.append(LenvmaxH)
       #vertical 
       AeventPtPZ.append(LptpZ)
       AeventEvelopeZ.append(LenvmaxZ)
    # change in array Horizonal
    AeventEvelopeH = np.asarray(AeventEvelopeH)
    AeventPtPH = np.asarray(AeventPtPH)
    
    # change in array vertical 
    AeventEvelopeZ = np.asarray(AeventEvelopeZ)
    AeventPtPZ = np.asarray(AeventPtPZ)
    
    np.savetxt(outfilePEnvH,AeventEvelopeH, fmt='%s')
    np.savetxt(outfilePtpH,AeventPtPH,fmt='%s')
    
    np.savetxt(outfilePEnvZ,AeventEvelopeZ, fmt='%s')
    np.savetxt(outfilePtpZ,AeventPtPZ,fmt='%s')
    
##########################################################################
#Phase 2 : Ratio and plot 
##########################################################################

LStations =np.delete(LStations, LStations.index(refStation))

Lcomponent = ['Z','H']
for Component in Lcomponent : 
    #ref array defintion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    outfilePtp = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_movingWindows_2Hz_0.2p_%s_%s.txt'%(refStation, Component)
    APtPref=np.loadtxt(outfilePtp)
    print APtPref.shape, 'Aptpref shape'
    fig,ax = plt.subplots()
    
    for Station in LStations: 
        outfilePtp = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_movingWindows_2Hz_0.2p_%s_%s.txt'%(Station, Component)
        APtP=np.loadtxt(outfilePtp)
        #if event not present in array UNIFORMISATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for i in xrange(APtPref.shape[0]):
             print APtPref[i][0:3], APtPref[i][0], 'tt'
             if i < len(APtP) :
                if  APtPref[i][0:3].all()<> APtP[i][0:3].all(): #[ligne][column]
                    #Julian day of ref is greater than the julian day of the non ref ==> add a day in the ref list 
                    if APtPref[i][0]>APtPref[i][0]:
                        print APtPref[i][0],APtP[i][0], 'l'
                        APtPref = np.insert(APtPref,i, APtP[i], axis=0)
                        APtPref[i][3]=np.nan
                    elif APtPref[i][0]<APtPref[i][0] :  #Julian day of r the other is lower than the julian day of the ref ==> add a day in the other list 
                        APtPref = np.insert(APtP,i, APtPref[i], axis=0)
                        APtP[i][3]=np.nan
                        #same day but maybe not same houror magnitude
                    else : 
                        #Hour of the ref day is greater than the hour of the non ref ==> add the line in the ref klist 
                        if APtPref[i][1]>APtPref[i][1]:
                            APtPref = np.insert(APtPref,i, APtP[i], axis=0)
                            APtPref[i][3]=np.nan
                        elif APtPref[i][1]<APtPref[i][1]:
                            APtPref = np.insert(APtP,i, APtPref[i], axis=0)
                            APtP[i][3]=np.nan
                        else : 
                            #HMagnitude of the ref day is greater than the hour of the non ref ==> add the line in the ref klist 
                            if APtPref[i][2]>APtPref[i][2]:
                                APtPref = np.insert(APtPref,i, APtP[i], axis=0)
                                APtPref[i][3]=np.nan
                            else: 
                                APtPref = np.insert(APtP,i, APtPref[i], axis=0)
                                APtP[i][3]=np.nan
             else :
                break
                    
        APtPcopy = np.delete(APtP,[0,1,2],axis=1)
        APtPrefcopy = np.delete(APtPref,[0,1,2],axis=1)
        Lratio = [APtPcopy[i]/APtPrefcopy[i] for i in xrange(APtPref.shape[0])]
        print len(Lratio), 'Lratio shape', APtPcopy, APtP, APtPref
        Lmean = np.nanmean(Lratio, axis =1)
        Lstd= np.nanstd(Lratio, axis =1)
        print len(Lmean)
        x= range(1.05,21,2)
        ax.plot(x,Lmean,color= Dcolor[Station], label ='S'+station )
        AMean.append(Lmean.tolist())
        
        Astd.append(Lsdt.tolist())
        
    plt.ylabel('Peak to Peak')
    plt.xlabel('Frequency (Hz)')       
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    outfile = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_movingWindows_refstation_%s_2Hz_0.2p_%s.eps'%(refStation, Component)
    plt.savefig(outfile,bbox_inches="tight")


##########################################################################
#Envelope 
##########################################################################

Lcomponent = ['Z','H']
for Component in Lcomponent : 
    #ref array defintion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    outfileEnv = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_movingWindows_2Hz_0.2p_%s_%s.txt'%(refStation, Component)
    AEnvref=np.loadtxt(outfileEnv)
    print AEnvref.shape, 'AEnvref shape'
    fig,ax = plt.subplots()
    
    for Station in LStations: 
        outfileEnv = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_movingWindows_2Hz_0.2p_%s_%s.txt'%(Station, Component)
        AEnv=np.loadtxt(outfileEnv)
        #if event not present in array UNIFORMISATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for i in xrange(AEnvref.shape[0]):
             print AEnvref[i][0:3], AEnvref[i][0], 'tt'
             if i < len(AEnv) :
                if  AEnvref[i][0:3].all()<> AEnv[i][0:3].all(): #[ligne][column]
                    #Julian day of ref is greater than the julian day of the non ref ==> add a day in the ref list 
                    if AEnvref[i][0]>AEnvref[i][0]:
                        print AEnvref[i][0],AEnv[i][0], 'l'
                        AEnvref = np.insert(AEnvref,i, AEnv[i], axis=0)
                        AEnvref[i][3]=np.nan
                    elif AEnvref[i][0]<AEnvref[i][0] :  #Julian day of r the other is lower than the julian day of the ref ==> add a day in the other list 
                        AEnvref = np.insert(AEnv,i, AEnvref[i], axis=0)
                        AEnv[i][3]=np.nan
                        #same day but maybe not same houror magnitude
                    else : 
                        #Hour of the ref day is greater than the hour of the non ref ==> add the line in the ref klist 
                        if AEnvref[i][1]>AEnvref[i][1]:
                            AEnvref = np.insert(AEnvref,i, AEnv[i], axis=0)
                            AEnvref[i][3]=np.nan
                        elif AEnvref[i][1]<AEnvref[i][1]:
                            AEnvref = np.insert(AEnv,i, AEnvref[i], axis=0)
                            AEnv[i][3]=np.nan
                        else : 
                            #HMagnitude of the ref day is greater than the hour of the non ref ==> add the line in the ref klist 
                            if AEnvref[i][2]>AEnvref[i][2]:
                                AEnvref = np.insert(AEnvref,i, AEnv[i], axis=0)
                                AEnvref[i][3]=np.nan
                            else: 
                                AEnvref = np.insert(AEnv,i, AEnvref[i], axis=0)
                                AEnv[i][3]=np.nan
             else :
                break
                    
        AEnvcopy = np.delete(AEnv,[0,1,2],axis=1)
        AEnvrefcopy = np.delete(AEnvref,[0,1,2],axis=1)
        Lratio = [AEnvcopy[i]/AEnvrefcopy[i] for i in xrange(AEnvref.shape[0])]
        print len(Lratio), 'Lratio shape', AEnvcopy, AEnv, AEnvref
        Lmean = np.nanmean(Lratio, axis =1)
        Lstd= np.nanstd(Lratio, axis =1)
        print len(Lmean)
        x= range(1.05,21,2)
        ax.plot(x,Lmean,color= Dcolor[Station], label ='S'+station )
        AMean.append(Lmean.tolist())
        
        Astd.append(Lsdt.tolist())
        
    plt.ylabel('Peak to Peak')
    plt.xlabel('Frequency (Hz)')       
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    outfile = '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/Envelope_movingWindows_refstation_%s_2Hz_0.2p_%s.eps'%(refStation, Component)
    plt.savefig(outfile,bbox_inches="tight")            
        
