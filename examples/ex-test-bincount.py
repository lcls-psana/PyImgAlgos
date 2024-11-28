#!/usr/bin/env python

import os
import sys
from time import time

import numpy as np
import pyimgalgos.GlobalUtils as gu
from pyimgalgos.HBins import HBins

shape = (1024,1024)

t0_sec = time()
arr = 100 + 10*np.random.standard_normal(size=shape).flatten()
print('\n np.random.standard_normal time %.6f sec' % (time()-t0_sec)) # 0.052328 sec
gu.print_ndarr(arr, 'np.random.standard_normal')


t0_sec = time()
arr_exp = 100*np.random.standard_exponential(size=shape).flatten()
arr_i16 = np.require(arr_exp, dtype=np.int16)
arr_int = np.require(arr_exp, dtype=np.int32)
arr_f32 = np.require(arr_exp, dtype=np.float32)
arr_f64 = np.require(arr_exp, dtype=np.float64)
print('\n np.random.standard_exponential time %.6f sec' % (time()-t0_sec)) # 0.052328 sec
gu.print_ndarr(arr_exp, 'np.random.standard_exponential')
gu.print_ndarr(arr_i16, '...i16')
gu.print_ndarr(arr_int, '...i  ')
gu.print_ndarr(arr_f32, '...f32')
gu.print_ndarr(arr_f64, '...f64')


hb = HBins((50,150), nbins=200)
t0_sec = time()
inds = hb.bin_indexes(arr) #, edgemode=0)
print('\n hb.bin_indexes time %.6f sec' % (time()-t0_sec)) # 0.010366 sec
gu.print_ndarr(inds, 'hb.bin_indexes')

#================


t0_sec = time()
bc = np.bincount(inds, weights=None)
print('\nnp.bincount inds time %.6f sec' % (time()-t0_sec)) # 0.002157 sec
gu.print_ndarr(bc, 'bincount:', last=100)



t0_sec = time()
bc = np.bincount(inds, weights=arr_f64)
print('\nnp.bincount f64 time %.6f sec' % (time()-t0_sec)) # 0.002157 sec
#gu.print_ndarr(bc, 'bincount:', last=100)


t0_sec = time()
bc = np.bincount(inds, weights=arr_f32)
print('np.bincount f32 time %.6f sec' % (time()-t0_sec)) # 0.002157 sec
#gu.print_ndarr(bc, 'bincount:', last=100)


t0_sec = time()
bc = np.bincount(inds, weights=arr_int)
print('np.bincount int time %.6f sec' % (time()-t0_sec)) # 0.002157 sec
#gu.print_ndarr(bc, 'bincount:', last=100)


t0_sec = time()
bc = np.bincount(inds, weights=arr_i16)
print('np.bincount i16 time %.6f sec' % (time()-t0_sec)) # 0.002157 sec
#gu.print_ndarr(bc, 'bincount:', last=100)

#------------------------------
