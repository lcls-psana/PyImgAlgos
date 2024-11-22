####!/usr/bin/env python

#--------------------------------
"""`pyimgalgos.MedianFilter.py` - contains method median_filter_ndarr and associated code for self test.

Usage::

    # Import
    # ======
    from pyimgalgos.MedianFilter import median_filter_ndarr

    #nda = ... should be defined somehow, for example nda = det.calib(evt)
    rank=3
    nda_med = median_filter_ndarr(nda, rank)

This software was developed for the SIT project.
If you use all or part of it, please give an appropriate acknowledgment.

Command to test python pyimgalgos/src/MedianFilter.py 3

Also see for details:
http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.ndimage.filters.median_filter.html#scipy.ndimage.filters.median_filter

Revision: $Revision$

@version $Id$

@author Mikhail S. Dubrovin

"""
from __future__ import print_function
from __future__ import division

#--------------------------------
__version__ = "$Revision$"
#--------------------------------

import sys
import numpy as np
from scipy.ndimage.filters import median_filter

#------------------------------

def reshape_nda_to_3d(arr) :
    """Reshape np.array to 3-d
    """
    sh = arr.shape
    if len(sh)<4 : return arr
    arr.shape = (arr.size//sh[-1]//sh[-2], sh[-2], sh[-1])
    return arr

#------------------------------

def footprint_ring(rank=3) :
    """Ring footprint for scipy.ndimage.filters.median_filter
       rank : int - radial parameter, pixels within radius less or equal rank will be used to evaluate the median.
    """
    r1 = rank+1
    shape_q = (r1,r1)
    grid = np.indices(shape_q)
    rank2 = rank**2
    fpq = np.array([r**2+c**2<=rank2 for r,c in zip(grid[0].flatten(),grid[1].flatten())], dtype=np.bool)
    fpq.shape = shape_q
    fph = np.hstack([fpq[:,::-1], fpq[:,1::]])
    fpr = np.vstack([fph[::-1,:], fph[1::,:]])
    #print fpq
    #print fpq[:,1::]
    #print '\n', fph
    #print '\n', fpr
    return fpr

#------------------------------

def median_filter_ndarr(nda_in, rank=3) :
    """returns 2-d or 3-d array with number of merged photons per pixel.

    Parameters
    - nda_in : numpy.array - n-dimensional numpy array
    - rank   : int - radial parameter, pixels within radius less or equal rank will be used to evaluate the median.

    use scipy.ndimage.filters.median_filter(input, size=None, footprint=None, output=None, mode='reflect', cval=0.0, origin=0)
    http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.ndimage.filters.median_filter.html#scipy.ndimage.filters.median_filter
    """
    shape_in = nda_in.shape

    ndim = len(shape_in)
    if ndim == 1 :
        msg = 'WARNING: pyimgalgos.median_filter_ndarr got 1-d array, while it works with ndim>1 imaging arrays - return None.'
        print(msg)
        #raise ValueError(msg)
        return None

    fp = footprint_ring(rank)
    if ndim == 2 : return median_filter(nda_in, size=None, footprint=fp, mode='mirror', origin=0)

    nda = nda_in if ndim==3 else reshape_nda_to_3d(nda_in)

    nda_out = np.array([median_filter(nda_in[s,:,:], footprint=fp, mode='mirror', origin=0)\
                        for s in range(nda.shape[0])], dtype=nda_in.dtype)
    nda_in.shape  = shape_in
    nda_out.shape = shape_in
    return nda_out

#------------------------------
#------------------------------
#----------- TEST -------------
#------------------------------
#------------------------------

def random_standard_ndarr(shape=(185,388), mu=50, sigma=10, dtype=np.float32) :
    """Returns n-d array of specified shape with random intensities generated for Gaussian parameters.
    """
    return (mu + sigma*np.random.standard_normal(shape)).astype(dtype,copy=False)

#------------------------------

def slope_2darr(shape=(185,388), axis=0, dtype=np.float32) :
    """Returns n-d array monotonicly raising along axis=0 or 1"""
    imax,jmax = shape
    if axis==0 :
        arr1d = list(range(imax))
        return np.array([arr1d for j in range(jmax)])
    else :
        arr1d = np.ones(imax)
        return np.array([arr1d*j for j in range(jmax)])

#------------------------------

def example01() :
    print("""example01 footprint_ring()""")
    fpr = footprint_ring(rank=50)
    if True :
      import pyimgalgos.GlobalGraphics as gg
      gg.plotImageLarge(fpr, amp_range=None, title='footprint')
      gg.show()

#------------------------------

def example02() :
    print("""example02 random_standard_ndarr(shape=(40,40))""")

    from time import time
    #nda = slope_2darr(shape=(5,5), axis=0, dtype=np.float32)
    nda = random_standard_ndarr(shape=(40,40), mu=50, sigma=10, dtype=np.float32)
    nda[10:30,10:30] = np.zeros((20,20))
    rank=3
    t0_sec = time()
    nda_med = median_filter_ndarr(nda, rank)
    print('median_filter_ndarr consumes time %.3f sec at rank=%d and array shape=%s' % (time()-t0_sec, rank, str(nda.shape)))

    if True :
      import pyimgalgos.GlobalGraphics as gg
      gg.plotImageLarge(nda, amp_range=None, title='nda_in')
      gg.plotImageLarge(nda_med, amp_range=None, title='median')
      gg.show()

#------------------------------

def example03() :
    print("""example03 random_standard_ndarr(shape=(32,185,388))""")

    from time import time
    nda = random_standard_ndarr(shape=(32,185,388), mu=50, sigma=10, dtype=np.float32)
    rank=3
    t0_sec = time()
    nda_med = median_filter_ndarr(nda, rank)
    print('median_filter_ndarr consumes time %.3f sec at rank=%d and array shape=%s' % (time()-t0_sec, rank, str(nda.shape)))

    if True :
      seg=1
      import pyimgalgos.GlobalGraphics as gg
      gg.plotImageLarge(nda[seg,:,:], amp_range=None, title='nda_in')
      gg.plotImageLarge(nda_med[seg,:,:], amp_range=None, title='median')
      gg.show()

#------------------------------

def usage() : return 'Use command: python %s <test-number [1-3]>' % sys.argv[0]

def main() :
    print('\n%s\n' % usage())
    if len(sys.argv)!= 2  : example01()
    elif sys.argv[1]=='1' : example01()
    elif sys.argv[1]=='2' : example02()
    elif sys.argv[1]=='3' : example03()
    else                  : sys.exit ('Test number parameter is not recognized.\n%s' % usage())

#------------------------------

if __name__ == "__main__" :
    main()
    sys.exit('End of test')

#------------------------------
