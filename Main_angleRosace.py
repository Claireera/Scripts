# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 21:40:46 2016

@author: claire
"""

"""main Angle"""
from Plot_results import *
from Events_Selections import *
station ='3'
File = '/home/claire/PHD/Working/Data/Event_stations_Caracteristics.txt'
MlMin=3
MlMax=9
DistMin=0
DistMax =200 
output = False
Aevent, AeventShape =SelectEventMLDist(File, station, MlMin, MlMax, DistMin, DistMax, output)
Langles = [Aevent[i][9] for i in xrange(Aevent.shape[0])]  #column 9 and all ligne
Lml =[Aevent[i][5] for i in xrange(Aevent.shape[0])]
print Aevent.shape
print len(Langles)

Anglediagram(Langles, '/home/claire/PHD/Working/Results/Rose_diagrams_events2015_2016_entiers4.eps')


#Histograms(Lml, len(xrange(int(min(Lml)),int(max(Lml))+1,1)),(min(Lml),max(Lml)+1),'magnitude', 'events', '/home/claire/PHD/Working/Results/Magnitude_histogram_2015_2016.eps')
