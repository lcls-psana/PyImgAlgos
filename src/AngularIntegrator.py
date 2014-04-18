#!/usr/bin/env python
#------------------------------
"""Evaluate radial histogram for angle-integrated 2-d image

This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@version $Id: 2008-09-28$

@author Mikhail S. Dubrovin
"""

#------------------------------
#  Module's version from CVS --
#------------------------------
__version__ = "$Revision$"
# $Source$

#--------------------------------

from time import time
import math
import numpy as np
#import scipy as sp
#import scipy.ndimage

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

#------------------------------

def valueToIndexProtected(V, VRange) :
    """Input: V - numpy array of values,
              VRange - Vmin, Vmax, contains the binning parameters Vmin, Vmax, and Nbins
       Output: Array of indexes from 0 to Nbins (Nbins+1 index) with shape of V.
       The last index Nbins is used for overflow and underflow.
    """
    Vmin, Vmax, Nbins = VRange
    Nbins1 = int(Nbins)-1
    factor = float(Nbins) / float(Vmax-Vmin)
    indarr = np.int32( factor * (V-Vmin) )
    #return np.select([V<Vmin,V>Vmax], [0,Nbins-1], default=indarr)
    return np.select([V==Vmax, indarr<0, indarr>Nbins1], [Nbins1, 0, 0], default=indarr)

#------------------------------

def divideArraysSafely(num, den) :
    """Per evement divides numpy arrays result = num/den. Protected for 0 values. Arrays should have the same size."""
    if num.shape != den.shape :
        print 'divideArraysSafely: non-equal array shapes for numerator and denumerator: ', num.shape, den.shape
    num_corr =  np.select([den<1], [0], default=num)    
    den_corr =  np.select([den<1], [1], default=den)    
    return num_corr/den_corr

#------------------------------

class AngularIntegrator :
    """Angular integration of a 2D numpy array"""
    
    def __init__(self):
        self.center  = None
        self.rrange  = None
        self.mask    = None
        self.imshape = None
        self.rpixind = None
        self.rbinind = None
        self.binstat = None

        
    def setParameters(self, xsize, ysize, xc=None, yc=None, rmin=None, rmax=None, nbins=None, mask=None):
        """Sets image, radial histogram parameters, and mask if proper normalization on number of actual pixels is desired
        """   
        # agregate input parameters
        self.mask    = mask
        self.imshape = (ysize, xsize)
        self.xsize   = xsize
        self.ysize   = ysize

        self.xc    = xsize/2                            if xc    is None else xc
        self.yc    = ysize/2                            if yc    is None else yc
        self.rmin  = 0                                  if rmin  is None else rmin
        self.rmax  = math.sqrt(self.xc**2 + self.yc**2) if rmax  is None else rmax
        self.nbins = int((self.rmax - self.rmin)/2)     if nbins is None else nbins

        self.center = (self.xc, self.yc)
        self.rrange = (self.rmin, self.rmax, self.nbins)

        # produce all index arrays
        x = np.arange(self.xsize) - self.xc
        y = np.arange(self.ysize) - self.yc
        xgrid, ygrid = np.meshgrid(x,y)
 
        half_bin = 0.5*(self.rmax-self.rmin)/self.nbins
        self.bincent = np.linspace(self.rmin+half_bin, self.rmax-half_bin, self.nbins)

        self.rpixind  = np.sqrt(xgrid**2 + ygrid**2)
        self.rbinind  = valueToIndexProtected(self.rpixind, self.rrange)

        if not self.mask is None : 
            self.binstat = np.bincount(self.rbinind.flatten(), weights=self.mask.flatten(), minlength=self.rrange[2])
            #print 'self.binstat   = ', self.binstat

        #print 'self.rpixind.shape =', self.rpixind.shape
        #print 'self.rbinind.shape =', self.rbinind.shape
        #print 'self.rbincent =',       self.bincent
        #print 'self.rbincent.shape =', self.bincent.shape


    def getRadialHistogramArrays(self, image):
        """Fills radial histogram with image intensities and do normalization on actual pixel number if the mask is provided
           and returns two arrays with radial bin centers and integrated (normalized) intensities. 
        """   
        bin_integral = np.bincount(self.rbinind.flatten(), weights=image.flatten(), minlength=self.rrange[2])
        #print 'bin_integral = ', bin_integral
        #print 'bin_integral.shape = ', bin_integral.shape
        if self.mask is None : return self.bincent, bin_integral

        return self.bincent, divideArraysSafely(bin_integral, self.binstat)


    def getRPixelIndexMap(self) :
        return self.rpixind


    def getRBinIndexMap(self) :
        return self.rbinind

#------------------------------
#------------------------------
#-------- FOR TEST ONLY -------
#------------------------------
#------------------------------

def drawImage(arr, img_range=None, zrange=None, figsize=(10,10)) :    # range = (left, right, low, high), zrange=(zmin,zmax)
    fig = plt.figure(figsize=figsize, dpi=80, facecolor='w', edgecolor='w', frameon=True)
    fig.subplots_adjust(left=0.05, bottom=0.03, right=0.98, top=0.98, wspace=0.2, hspace=0.1)
    figAxes = fig.add_subplot(111)
    #figAxes = fig.add_axes([0.15, 0.06, 0.78, 0.21])
    imAxes = figAxes.imshow(arr, origin='upper', interpolation='nearest', aspect='auto', extent=img_range)
    if zrange != None : imAxes.set_clim(zrange[0],zrange[1])
    colbar = fig.colorbar(imAxes, pad=0.03, fraction=0.04, shrink=1.0, aspect=40, orientation='horizontal')

#------------------------------

def drawGraph(x,y) : 
    fig = plt.figure(figsize=(6,4), dpi=80, facecolor='w', edgecolor='w', frameon=True)
    #fig.subplots_adjust(left=0.05, bottom=0.03, right=0.98, top=0.98, wspace=0.2, hspace=0.1)
    #figAxes = fig.add_subplot(111)
    ax = fig.add_axes([0.15, 0.10, 0.78, 0.86])
    ax.plot(x,y,'b-')

#------------------------------

def MakeImage(shape=(1024,1024)) :
    # Create test image - a sinc function, centered in the middle of
    # the image
    # Integrating in phi about center, the integral will become sin(x)
    print "Creating test image",

    xsize, ysize = shape
    ratio = float(ysize)/float(xsize)
    print 'ratio = ', ratio
    xmin, xmax = -4, 6
    ymin, ymax = -7*ratio, 3*ratio

    print '\nxmin, xmax, xsize = ', xmin, xmax, xsize
    print '\nymin, ymax, ysize = ', ymin, ymax, ysize

    xarr = np.linspace(xmin, xmax, xsize)
    yarr = np.linspace(ymin, ymax, ysize)
    xgrid, ygrid = np.meshgrid(xarr, yarr)
    rgrid = np.sqrt(xgrid**2 + ygrid**2)
    image = np.abs(np.sinc(rgrid))
    return image

#------------------------------

def mainTest() :

    image = MakeImage(shape=(2000,1000))
    mask  = np.ones_like(image, dtype=np.int)

    ysize, xsize = image.shape

    t0_sec = time()
    angint = AngularIntegrator()
    angint.setParameters(xsize, ysize, xc=xsize*0.4, yc=ysize*0.7, rmin=0, rmax=1000, nbins=1000, mask=mask)
    #angint.setParameters(xsize, ysize, rmin=1, rmax=1200, nbins=200, mask=mask)
    print "Time consumed for indexing = %.3f sec" % (time()-t0_sec)

    drawImage(angint.getRBinIndexMap(), figsize=(10,6))

    drawImage(image, figsize=(10,6))

    t0_sec = time()
    bincent, integral = angint.getRadialHistogramArrays(image)
    print "Time consumed for intensity binning = %.3f sec" % (time()-t0_sec)

    #plt.plot(bincent, integral)
    drawGraph(bincent, integral)

    plt.show()
 
#------------------------------

if __name__ == "__main__" :

    mainTest()

#------------------------------
