#!/usr/bin/env python

#------------------------------
from __future__ import print_function
import sys
from time import time
import numpy as np
from pyimgalgos.GlobalUtils import print_ndarr
#------------------------------

shape = (1000,1000)
mean, rms = 100,10 
arr  = np.array(mean + rms * np.random.standard_normal(size=shape), dtype=np.int16)

sta  = np.zeros(shape, dtype=np.int16)
ones = np.ones(shape, dtype=np.int16)

print_ndarr(arr, 'arr', first=0, last=20)

#------------------------------

t0_sec = time()
conds = arr < (mean-0.5*rms)
print('Time = %9.6f sec' % (time()-t0_sec))
print_ndarr(conds, 'conds', first=0, last=20)
# Time =  0.001355 sec

#------------------------------

t0_sec = time()
sta[conds] += 1
print('Time = %9.6f sec' % (time()-t0_sec))
print_ndarr(sta, 'statistics1', first=0, last=20)
# Time =  0.010388 sec

#------------------------------

t0_sec = time()
sta[conds] += ones[conds]
print('Time = %9.6f sec' % (time()-t0_sec))
print_ndarr(sta, 'statistics2', first=0, last=20)
#Time =  0.015995 sec

#------------------------------

t0_sec = time()
sta += np.select((conds,), (ones,), 0)
print('Time = %9.6f sec' % (time()-t0_sec))
print_ndarr(sta, 'statistics3', first=0, last=20)
#Time =  0.006738 sec

#------------------------------

t0_sec = time()
sta += np.select((conds,), (1,), 0)
print('Time = %9.6f sec' % (time()-t0_sec))
print_ndarr(sta, 'statistics4', first=0, last=20)
#Time =  0.008723 sec

#------------------------------
#------------------------------
sys.exit('END OF TEST')
#------------------------------
#------------------------------

