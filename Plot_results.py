# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 23:19:11 2016

@author: claire
"""
"""funtions :
Anglediagram :  from a List of angle return a plot diagrams of pourcentages of event between a angle thereshold 
Anglediagram_entier:from a List of angle return a plot diagrams of the number of event between a angle thereshold 
Histograms : histograms of a given array return plot, histogram and bin_edges 
phase : plot phase 
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.colorbar import ColorbarBase
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from operator import itemgetter
import math

    
def Anglediagram(Langles, outputfile):
    """ return a rose diagram of the proportion of events between  a given azimuth thershold [0,10[
   *  input : 
    -Langle :  List of event azimuthn in degrees [-180,180] : [-180;0] East, [0;180] Ouest 
    - outputfile :  directory and file name 
    * output : 
        rose diagram with pourcentages of events between given azimuth save in outputfile
    * exemple :
     Langles = [1,2,56,-158,160,-12]
     Anglediagram(Langles, '/home/claire/PHD/Working/Results/Rose_diagrams_events2015_2016.eps') """
    #build the rosace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    l = len(Langles)
    angle = math.radians(10.)
    patches = math.radians(360.)/angle
    theta = np.arange(0,math.radians(360.),angle)
    count = [0]*int(patches)
    # angle classification and counting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i, angi in enumerate(Langles):
        if angi<0:
            #angle is given between -180;180 
            ang=360+angi
        else:
            ang=angi
        # azimuth conversion in radian ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        item=math.radians(ang)
        #counting per thereshold~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        temp = int((item - item%angle)/angle)
        count[temp] += 1
        
    #width of the cololum as a function of the thershiold width~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    width = angle * np.ones(patches)
    # force square figure and square axes looks better for polar, IMO
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    countpercent =[a*100/l for a in count] 
    rmax = max(countpercent) + 1
    #set axis label, grid ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ax.set_rlim(0,rmax)
    ax.set_theta_offset(np.pi/2)
    ax.set_thetagrids(np.arange(0,360,10))
    ax.set_theta_direction(-1)
    
    # project strike distribution as histogram bars
   
    bars = ax.bar(theta,countpercent, width=width)
    r_values = []
    colors = []
    for r,bar in zip(countpercent, bars):
        r_values.append(r/float(max(countpercent)))
        #choose the last color input wich is proportional to the number of event considered~~~~~~~~~~~~~~~~
        colors.append(cm.Blues(r_values[-1]))
        bar.set_facecolor(colors[-1])
        #bar.set_edgecolor('grey')
        #bar.set_alpha(0.5)
    
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    colorlist = []
    r_values.sort()
    values = []
    for val in r_values:
        if val not in values:
            values.append(val*float(max(count)))
    
        color = cm.Blues(val)
        if color not in colorlist:
            colorlist.append(color)
    
#    cpt = ListedColormap(colorlist)
#    
    plt.show()
    plt.savefig(outputfile,bbox_inches="tight")
    return

def Anglediagram_entier(Langles, outputfile):
    """ return a rose diagram of the proportion of events between  a given azimuth thershold [0,10[
   *  input : 
    -Langle :  List of event azimuthn in degrees [-180,180] : [-180;0] East, [0;180] Ouest 
    - outputfile :  directory and file name 
    * output : 
        rose diagram with pourcentages of events between given azimuth save in outputfile
    * exemple :
     Langles = [1,2,56,-158,160,-12]
     Anglediagram(Langles, '/home/claire/PHD/Working/Results/Rose_diagrams_events2015_2016.eps') """
    #build the rosace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    l = len(Langles)
    angle = math.radians(10.)
    patches = math.radians(360.)/angle
    theta = np.arange(0,math.radians(360.),angle)
    count = [0]*int(patches)
    # angle classification and counting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i, angi in enumerate(Langles):
        if angi<0:
            #angle is given between -180;180 
            ang=360+angi
        else:
            ang=angi
        # azimuth conversion in radian ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        item=math.radians(ang)
        #counting per thereshold~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        temp = int((item - item%angle)/angle)
        count[temp] += 1
        
    #width of the cololum as a function of the thershiold width~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    width = angle * np.ones(patches)
    # force square figure and square axes looks better for polar, IMO
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    
    rmax = max(count) + 1
    #set axis label, grid ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ax.set_rlim(0,rmax)
    ax.set_theta_offset(np.pi/2)
    ax.set_thetagrids(np.arange(0,360,10))
    ax.set_theta_direction(-1)
    
    # project strike distribution as histogram bars
   
    bars = ax.bar(theta,count, width=width)
    r_values = []
    colors = []
    for r,bar in zip(count, bars):
        r_values.append(r/float(max(count)))
        #choose the last color input wich is proportional to the number of event considered~~~~~~~~~~~~~~~~
        colors.append(cm.Blues(r_values[-1]))
        bar.set_facecolor(colors[-1])
        #bar.set_edgecolor('grey')
        #bar.set_alpha(0.5)
    
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    colorlist = []
    r_values.sort()
    values = []
    for val in r_values:
        if val not in values:
            values.append(val*float(max(count)))
    
        color = cm.Blues(val)
        if color not in colorlist:
            colorlist.append(color)
    
#    cpt = ListedColormap(colorlist)
#    
    plt.show()
    plt.savefig(outputfile,bbox_inches="tight")
    return


def Histograms(L, binNb,Range, Xlabel, Ylabel, outfile):
    """ L type array, array of data to count
    * input :
        - L :  type np.array ; Data array to plot 
        - binNb number of bins 
        - Range : L.min and L.max() would be taken if no data is given ex None or [0.2,20]
        - Xlabel : type str; name of the X axes
        - Ylabel : type str; name of the Y axes
        - outfile : type str ; directory and name of the output file 
    * output :
    - hist :  type array ; values of histograms 
    - bin_edges : type arrays; bin edges 
    * exemple :
    L= np.array([1,2,6,7,5,5,5,7,9])
    Histograms(L,len(xrange(min(L),max(L)+1,1)),(min(L),max(L)+1), 'magnitude', 'events', '/home/claire/PHD/Working/Results/Magnitude_histogram_2015_2016.eps')
    """
    Hist, bin_edges=np.histogram(L, bins=binNb, range=Range, normed=False, weights=None, density=None)
    
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    B= [(bin_edges[i]+bin_edges[i+1])/2 for i in xrange(len(bin_edges)-1)]
    np.append(Hist,0)
    Hist = [float(a)/np.sum(Hist) for a in Hist]
    ax.bar(B,Hist,align='center', width = 1, color= '#E6E6E6', alpha =0.5)
    plt.savefig(outfile,bbox_inches="tight")
    return Hist, bin_edges
    
    
def ParticulemotionDiagrams(trZ, trH, step, plot, outfile):
    """plot velocity particule motion by plotting the vertical as the Horizontal component of the motion
    * input
        - trZ: type trace, vertical trace
        - trH: type trace, horizontal trace 
        - step: plot motion every step in number of points
    * output
    * exemple : 
    outfile = '/home/claire/PHD/Working/Results/Particule_motion/Plots_Particules_Motion/ParticuleMotion_Year_Jjul_Second.eps'
     ParticulemotionDiagrams(trZ, trH, step, plot, outfile)
    """
    H=[]
    Z=[]
    #record trz and trh valies every step (10points)
    for j in xrange(0, len(trZ), step):
        H.append(trH.data[j])
        Z.append(trZ.data[j])
    plt.subplot(111)
    plt.plot(H,Z, 'k')
    plt.xlabel('$Horizontal \; m.s^-1$')
    plt.ylabel('$Vertical \; m.s^-1$')
    plt.savefig(outfile,bbox_inches="tight")
    return 
    
def MeanSTDHVSR(AHVSR,Afr, outfile):
    Mean = np.mean(AHVSR,axis=0)
    fr = np.mean(Afr, axis=0)
    print Mean
    STD = np.std(AHVSR,axis=0)
    STD_u = Mean + STD
    STD_l = Mean - STD
    plt.figure()
    plt.plot(fr,Mean,'k', label = 'Mean')
    plt.plot(fr,STD_u,'#787878', linestyle = '--',label = 'Std + ')
    plt.plot(fr,STD_l,'#787878',linestyle = '--', label = 'Std - ')
    plt.legend()
    plt.xlim(0.01,20)
    plt.ylim(-1,10)
    plt.ylabel("HVSR")
    plt.xlabel("Frequency Hz")
    plt.savefig(outfile)
    return 


def PlotSpectrogram(FrArray, HVSRArray, JjulArray, outfile,colormap="jet"):
    plt.figure()
    FrArray = [ FrArray[:][i] for i in xrange(0, len(FrArray),5)]
    FrArray=FrArray[:][0:len(FrArray)/2]
    print len(FrArray)
    print HVSRArray.shape
    HVSRArray = np.array([HVSRArray[:][i] for i in xrange(0,HVSRArray.shape[0],5)])
    HVSRArray = HVSRArray[:][0:HVSRArray.shape[0]/2]
    print HVSRArray.shape
    X, Y = np.meshgrid(JjulArray, FrArray)
    plt.pcolor(X, Y,HVSRArray, cmap=colormap, vmin=0, vmax= 10)
    plt.colorbar()

    plt.xlabel("time (s)")
    plt.ylabel("frequency (Hz)")
    plt.ylim([0, 20])
   
    plt.savefig(outfile, bbox_inches="tight")
   
    return 
    
def PeaktoPeakPlot(Lstation,Fmin,Fmax,outfile)  :
    """ return, plot of the ratio of peak to peak for a given reference station
    *input:
       - Lstation, tyep list; list of the station to consider
       -Fmin,type float; min frequency used in the passband filter
       - Fmax,type float; max frequency used in the passband filter
        outfile type str, name and folder of the output plot
    
    *output:
        plot save in output file
    *exemple :
        PeaktoPeakPlot(Lstation, refstation, outfile)
         PeaktoPeakPlot(['1','2','3','4','5','6'], '7', '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_PBafterQsum_Lptpref_%s.txt'%'7')
    """
    plt.figure()
    Dcolor={'1':'#2ba70f','2':'#3342ff','3':'#ff3333','4':'#ffac33','5':'#D7DF01','6':'#4B0082','7':'#FF00FF'}
    Station=0
    for Station in Lstation: 
        Lptp = np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_Lptp%s%s.txt'%(Station,str(Fmin),str(Fmax)))
        Jjul = np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUl%s%s.txt'%(Station,str(Fmin),str(Fmax)))
        print Lptp, Jjul
#        Lptp = np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_Lptp.txt'%Station)
#        Jjul = np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUl.txt'%Station) 
        plt.plot(Jjul, Lptp, color = Dcolor[Station], marker = '.', lw=0.2, label = 'S' +Station )
        plt.ylabel('Peak to Peak')
        plt.xlabel('Julian day')
        plt.ylim(ymax=0.0001)
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.savefig(outfile,bbox_inches="tight")
    return 
 

def RatioPeaktoPeakPlot(Lstation,Fmin,Fmax,refstation, outfile):
    """ return, plot of the ratio of peak to peak for a given reference station
    *input:
        - Lstation, tyep list; list of the station to consider
        - Fmin,type float; min frequency used in the passband filter
        - Fmax,type float; max frequency used in the passband filter
        - refstation, type str,; reference sttion name 
        - outfile type str, name and folder of the output plot
    *output:
        plot save in output file
    *exemple :
         RatioPeaktoPeakPlot(['1','2','3','4','5','6'], '7', '/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_PBafterQsum_Lptpref_%s.txt'%'7')
    """
    
    Dcolor={'1':'#2ba70f','2':'#3342ff','3':'#ff3333','4':'#ffac33','5':'#D7DF01','6':'#4B0082','7':'#FF00FF'}
    Station = refstation
    #Lptpref = np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_Lptp%s%st.txt'%(Station,str(Fmin),str(Fmax)))
    JUlLPtPref = np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUlMlPytp%s%st.txt'%(Station,str(Fmin),str(Fmax)))
#    Lptpref = np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_Lptp.txt'%refstation)
#    JUlref = np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUl.txt'%refstation) 
#    L = JUlLPtPref  
#    print L, 'll'
    JUlLPtPref=sorted(JUlLPtPref, key=itemgetter(0,1,2))
#    print len(L)
#    JUlref = [L[i][0] for i in xrange(len(L))]
#    print JUlref, 'julref'
#    Lptpref= [L[i][3] for i in xrange(len(L))]
#  
    plt.figure()
    for Station in Lstation:
        Lratio =[]
        LJul=[]
        JUlLPtP= np.loadtxt('/home/claire/PHD/Working/Results/Peak_to_Peak/PeaktoPeak_File/PeaktoPeak_%s_PBafterQsum_LjJUlMlPytp%s%s.txt'%(Station,str(Fmin),str(Fmax)))
 
#        L1 = JUlLPtP
#        print L1, 'll'
        JUlLPtP=sorted(JUlLPtP, key=itemgetter(0,1,2))
#        print len(L1)
#        Jjul = [L1[i][0] for i in xrange(len(L1))]
#        print JUlref, 'julref'
#        Lptp= [L1[i][3] for i in xrange(len(L1))]
        i = 0
#        a = len(JUlref)
        JUlLPtP=np.array([JUlLPtP[i].tolist() for i in xrange(len(JUlLPtPref))])
        for i in xrange(len(JUlLPtPref)):
#            a=a-1
#            print 'not same event'  
#            if i<len(Jjul):
            print JUlLPtPref[i][0:3], JUlLPtPref[i][0], 'tt'
            if i < len(JUlLPtP):
                if  JUlLPtPref[i][0:3].all()<> JUlLPtP[i][0:3].all(): #[ligne][column]
                    #Julian day of ref is greater than the julian day of the non ref ==> add a day in the ref list 
                    if JUlLPtPref[i][0]>JUlLPtPref[i][0]:
                        print JUlLPtPref[i][0],JUlLPtP[i][0], 'l'
                        JUlLPtPref = np.insert(JUlLPtPref,i, JUlLPtP[i], axis=0)
                        JUlLPtPref[i][3]=np.nan
                    elif JUlLPtPref[i][0]<JUlLPtPref[i][0] :  #Julian day of r the other is lower than the julian day of the ref ==> add a day in the other list 
                        JUlLPtPref = np.insert(JUlLPtP,i, JUlLPtPref[i], axis=0)
                        JUlLPtP[i][3]=np.nan
                        #same day but maybe not same houror magnitude
                    else : 
                        #Hour of the ref day is greater than the hour of the non ref ==> add the line in the ref klist 
                        if JUlLPtPref[i][1]>JUlLPtPref[i][1]:
                            JUlLPtPref = np.insert(JUlLPtPref,i, JUlLPtP[i], axis=0)
                            JUlLPtPref[i][3]=np.nan
                        elif JUlLPtPref[i][1]<JUlLPtPref[i][1]:
                            JUlLPtPref = np.insert(JUlLPtP,i, JUlLPtPref[i], axis=0)
                            JUlLPtP[i][3]=np.nan
                        else : 
                            #HMagnitude of the ref day is greater than the hour of the non ref ==> add the line in the ref klist 
                            if JUlLPtPref[i][2]>JUlLPtPref[i][2]:
                                JUlLPtPref = np.insert(JUlLPtPref,i, JUlLPtP[i], axis=0)
                                JUlLPtPref[i][3]=np.nan
                            else: 
                                JUlLPtPref = np.insert(JUlLPtP,i, JUlLPtPref[i], axis=0)
                                JUlLPtP[i][3]=np.nan
                    
                    
#                    Jjul = np.insert(Jjul,i, JUlref[i], axis=0)
#                    print Jjul,i, JUlref[i]
#                    Lptp = np.insert(Lptp,i, np.nan,axis=0)
#                    a=a+1
#                    break
#                elif JUlref[i]>Jjul[i]:
#                    JUlref = np.insert(JUlref,i,Jjul[i], axis=0)
#                    print JUlref[i],i,Jjul[i],
#                    Lptpref = np.insert(Lptpref,i, np.nan,axis=0)
#                    print 
#                    break
#                    
            else :
                break
        LJul = np.array([JUlLPtPref[i].tolist() for i in xrange(len(JUlLPtPref))])[:,0]
        print LJul, 'lljul'
        LPtPref = np.array([JUlLPtPref[i].tolist() for i in xrange(len(JUlLPtPref))])[:,3]
        LPtP = np.array([JUlLPtP[i].tolist() for i in xrange(len(JUlLPtP))])[:,3]
        Lratio= [LPtP[i]/ LPtPref[i] for i in range(np.min([len( LPtP),len(LPtPref)]))]
        LJul = np.array([[i] for i in LJul[0:len(Lratio)]])
        Lratio = np.array([[i] for i in Lratio])
        print  LJul.shape, Lratio.shape
        AJulRatios=np.array(np.hstack((LJul,Lratio)))
        print AJulRatios[:,0],[str(AJulRatios[i][0]) for i in xrange(AJulRatios.shape[0])]
        
        plt.plot(AJulRatios[:,0],AJulRatios[:,1], color = Dcolor[Station], marker = '.' ,lw=0.1, label = 'S'+Station)
    plt.xticks(AJulRatios[:,0],[str(AJulRatios[i][0]) for i in xrange(AJulRatios.shape[0])],rotation='vertical')   
    plt.subplots_adjust(bottom=0.15) #spacing to prevent clipped byof ticks labels 
    plt.ylabel('Peak to Peak')
    plt.xlabel('Julian day')
    plt.ylim(ymax=10)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.savefig(outfile,bbox_inches="tight")
    return LJul,Lratio