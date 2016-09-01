# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 23:37:36 2016

@author: claire
"""

"""Main HVSR Mean between Fmin and Fmax, time variation """

import numpyt as np 

JJul = range(59,365)
Fmin = 
Fmax = 
fr = np.load()
Afr = np.arwhere((fr>Fmin)&(fr<Fmax))
MinA = np.min(Afr)
MaxA =np.max(Afr)
for Station in Lstation:
    AHVSR = np.load()
    AHVSR = AHVSR[MinA:MaxA][:] #crop the matrix 
    MeanA = np.nanmean(A, axis=1)

    print len(MeanA), len(Jjul)

    plt.plot(JJul,MeanA, label = 'S'+Station)
plt.xlabel('Julian Day')
plt.xlabel('HVSR Average between %s and %s'%(str(Fmin), str(Fmax)))
plt.savefig()