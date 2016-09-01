# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 13:47:51 2016

@author: claire
"""

#test h5py


#from mpi4py import MPI
import h5py
#Comm = MPI.COMM_WORLD
#
#
#rank = Comm.Get_rank()  # The process ID (integer 0-3 for 4-process run)
#size = Comm.Get_size()  # The number of processes
EQname = '15_12_01_25'
Component ='N'
f1 = h5py.File('/home/claire/PHD/Working/parallel_test.hdf5', 'w')
grpEQ = f1.create_group(EQname)
print 'grp'
dset = grpEQ.create_dataset('PtP_{}'.format(Component),(11,), dtype='i')
dset.attrs['PGA'] = 0
f1.close()

