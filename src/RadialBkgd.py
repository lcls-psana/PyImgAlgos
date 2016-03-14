#!/usr/bin/env python

#------------------------------
""":py:class:`RadialBkgd` - radial background subtraction for imaging detector n-d array data

Usage::

    # Import
    # ------
    from pyimgalgos.RadialBkgd import RadialBkgd

    # Initialization
    # --------------
    rb = RadialBkgd(xarr, yarr, mask=None, radedges=None, nradbins=100, phiedges=(0,360), nphibins=32)

    # Access methods
    # --------------
    orb          = rb.obj_radbins() # returns HBins object
    opb          = rb.obj_phibins() # returns HBins object
    rad          = rb.pixel_rad()
    irad         = rb.pixel_irad()
    phi          = rb.pixel_phi()
    iphi         = rb.pixel_iphi()
    iseq         = rb.pixel_iseq()
    npix_per_bin = rb.npixels_per_bin()
    int_per_bin  = rb.intensity_per_bin(nda)
    bkgd         = rb.bkgd_nda(nda)
    res          = rb.subtract_bkgd(nda)

    # Global methods
    # --------------
    result = divide_protected(num, den, vsub_zero=0)
    r, theta = cart2polar(x, y)
    x, y = polar2cart(r, theta)
    bin_values = bincount(map_bins, map_weights=None, length=None)

@see :py:class:`pyimgalgos.RadialBkgd`
See `Radial background <https://confluence.slac.stanford.edu/display/PSDMInternal/Radial+background+subtraction+algorithm>`_.

This software was developed for the SIT project.
If you use all or part of it, please give an appropriate acknowledgment.

Revision: $Revision$

@version $Id$

@author Mikhail S. Dubrovin

"""
#--------------------------------
__version__ = "$Revision$"
#--------------------------------

import math
import numpy as np
from pyimgalgos.HBins import HBins

import pyimgalgos.GlobalUtils as gu

#------------------------------

def divide_protected(num, den, vsub_zero=0) :
    """Returns result of devision of numpy arrays num/den with substitution of value vsub_zero for zero den elements.
    """
    pro_num = np.select((den!=0,), (num,), default=vsub_zero)
    pro_den = np.select((den!=0,), (den,), default=1)
    return pro_num / pro_den


def cart2polar(x, y) :
    """For numpy arrays x and y returns the numpy arrays of r and theta 
    """
    r = np.sqrt(x*x + y*y)
    theta = np.rad2deg(np.arctan2(y, x)) #[-180,180]
    #theta0 = np.select([theta<0, theta>=0],[theta+360,theta]) #[0,360]
    return r, theta


def polar2cart(r, theta) :
    """For numpy arryys r and theta returns the numpy arrays of x and y 
    """
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


def bincount(map_bins, map_weights=None, length=None) :
    """Wrapper for numpy.bincount with protection of weights and alattening numpy arrays
    """
    weights = None if map_weights is None else map_weights.flatten()
    return np.bincount(map_bins.flatten(), weights, length)

#------------------------------

class RadialBkgd() :
    def __init__(self, xarr, yarr, mask=None, radedges=None, nradbins=100, phiedges=(0,360), nphibins=32) :
        """Parameters
           - mask     - n-d array with mask
           - xarr     - n-d array with pixel x coordinates in any units
           - yarr     - n-d array with pixel y coordinates in the same units as xarr
           - radedges - radial bin edges for corrected region in the same units of xarr;
                        default=None - all radial range
           - nradbins - number of radial bins
           - phiedges - phi ange bin edges for corrected region.
                        default=(0,360)
                        Difference of the edge limits should not exceed +/-360 degree 
           - nphibins - number of angular bins
                        default=32 - bin size equal to 1 rhumb for default phiedges
        """
        self.mask = mask
        self.rad, self.phi = cart2polar(xarr, yarr)

        phi0 = min(phiedges[0], phiedges[-1])
        print 'Phi limits: ', phiedges[0], phiedges[-1]

        self.phi = np.select((self.phi<phi0, self.phi>=phi0), (self.phi+360.,self.phi))

        if len(self.rad.shape)>3 :
            self.rad = gu.reshape_to_3d(self.rad)
            self.phi = gu.reshape_to_3d(self.phi)

        self._set_rad_bins(radedges, nradbins)
        self._set_phi_bins(phiedges, nphibins)
        
        npbins = self.pb.nbins()
        nrbins = self.rb.nbins()
        ntbins = npbins*nrbins
        
        self.irad = self.rb.bin_indexes(self.rad, edgemode=1)        
        self.iphi = self.pb.bin_indexes(self.phi, edgemode=1)

        cond = np.logical_and(\
               np.logical_and(self.irad > -1, self.irad < nrbins),
               np.logical_and(self.iphi > -1, self.iphi < npbins)
               )
        if mask is not None : cond *= gu.reshape_to_3d(mask)

        self.iseq = np.select((cond,), (self.iphi*nrbins + self.irad,), ntbins).flatten()

        self.npix_per_bin = np.bincount(self.iseq, weights=None, minlength=None)


    def _set_rad_bins(self, radedges, nradbins) :
        rmin = math.floor(np.amin(self.rad)) if radedges is None else radedges[0]
        rmax = math.ceil (np.amax(self.rad)) if radedges is None else radedges[-1]
        if rmin<1 : rmin = 1
        self.rb = HBins((rmin, rmax), nradbins)


    def _set_phi_bins(self, phiedges, nphibins) :
        if phiedges[-1] > phiedges[0]+360\
        or phiedges[-1] < phiedges[0]-360:
            raise ValueError('Difference between angular edges should not exceed 360 degree;'\
                             ' phiedges: %.0f, %.0f' % (phiedges[0], phiedges[-1]))        
        self.pb = HBins(phiedges, nphibins)


    def print_attrs(self) :
        print self.pb.strrange(fmt='Phi bins:  min:%6.0f  max:%6.0f  nbins:%5d')
        print self.rb.strrange(fmt='Rad bins:  min:%6.0f  max:%6.0f  nbins:%5d')


    def print_ndars(self) :
        from Detector.GlobalUtils import print_ndarr
        print_ndarr(self.rad,'self.rad')
        print_ndarr(self.phi,'self.phi')


    def pixel_obj_radbins(self) :
        return self.rb


    def pixel_obj_phibins(self) :
        return self.pb


    def pixel_rad(self) :
        return self.rad


    def pixel_irad(self) :
        return self.irad


    def pixel_phi(self) :
        return self.phi


    def pixel_iphi(self) :
        return self.iphi


    def pixel_iseq(self) :
        return self.iseq


    def npixels_per_bin(self) :
        return self.npix_per_bin


    def intensity_per_bin(self, nda) :
        return np.bincount(self.iseq, weights=nda.flatten(), minlength=None)


    def average_per_bin(self, nda) :
        num = self.intensity_per_bin(nda)
        den = self.npixels_per_bin()
        return divide_protected(num, den, vsub_zero=0)


    def bkgd_nda(self, nda) :
        bin_bkgd = self.average_per_bin(nda)
        return np.array([bin_bkgd[i] for i in self.iseq])


    def subtract_bkgd(self, nda) :
        return nda.flatten() - self.bkgd_nda(nda)

#------------------------------

def test(ntest) :

    from time import time
    #from Detector.GlobalUtils import print_ndarr
    import pyimgalgos.GlobalGraphics as gg
    from PSCalib.GeometryAccess import GeometryAccess, img_from_pixel_arrays
    from PSCalib.NDArrIO import save_txt, load_txt

    #from pyimgalgos.NDArrGenerators import random_standard
    #arr = random_standard(shape=(40,60), mu=200, sigma=25)

    #import psana
    #ds  = psana.DataSource('exp=cxij4716:run=22')
    #det = psana.Detector('CxiDs2.0:Cspad.0', ds.env())

    dir       = '/reg/g/psdm/detector/alignment/cspad/calib-cxi-camera2-2016-02-05'
    fname_nda = '%s/nda-water-ring-cxij4716-r0022-e000001-CxiDs2-0-Cspad-0-ave.txt' % dir
    #fname_nda = '%s/nda-water-ring-cxij4716-r0022-e001000-CxiDs2-0-Cspad-0-ave.txt' % dir
    fname_geo = '%s/calib/CsPad::CalibV1/CxiDs2.0:Cspad.0/geometry/geo-cxi01516-2016-02-18-Ag-behenate-tuned.data' % dir

    # load n-d array with averaged water ring
    arr = load_txt(fname_nda)
    #print_ndarr(arr,'water ring')
    arr.shape = (32*185, 388)

    # retrieve geometry
    t0_sec = time()
    geo = GeometryAccess(fname_geo)
    iX, iY = geo.get_pixel_coord_indexes()
    X, Y, Z = geo.get_pixel_coords()
    mask = geo.get_pixel_mask(mbits=0377)
    print 'Time to retrieve geometry %.3f sec' % (time()-t0_sec)

    t0_sec = time()
    rbkg = RadialBkgd(X, Y, mask) # v0
    #rbkg = RadialBkgd(X, Y, mask, nradbins=500, nphibins=1) # v1
    #rbkg = RadialBkgd(X, Y, mask, nradbins=500, nphibins=8, phiedges=(-20, 240), radedges=(10000,80000)) # v2
    #rbkg = RadialBkgd(X, Y, mask, nradbins=3, nphibins=8, phiedges=(240, -20), radedges=(80000,10000)) # v3
    #rbkg = RadialBkgd(X, Y, mask, nradbins=3, nphibins=8, phiedges=(-20, 240), radedges=(10000,80000))

    print 'RadialBkgd initialization time %.3f sec' % (time()-t0_sec)

    #print 'npixels_per_bin:',   rbkg.npixels_per_bin()
    #print 'intensity_per_bin:', rbkg.intensity_per_bin(arr)
    #print 'average_per_bin:',   rbkg.average_per_bin(arr)

    t0_sec = time()
    nda = arr
    if ntest == 2 : nda = rbkg.pixel_rad()
    if ntest == 3 : nda = rbkg.pixel_phi()
    if ntest == 4 : nda = rbkg.pixel_irad() + 2
    if ntest == 5 : nda = rbkg.pixel_iphi() + 2
    if ntest == 6 : nda = rbkg.pixel_iseq() + 2
    if ntest == 7 : nda = mask
    if ntest == 8 : nda = rbkg.bkgd_nda(nda)
    if ntest == 9 : nda = rbkg.subtract_bkgd(nda) * mask.flatten() 
    print 'Get n-d array time %.3f sec' % (time()-t0_sec)

    img = img_from_pixel_arrays(iX, iY, nda)

    da, ds = None, None

    if ntest in (2,3,4,5,6,7) :
      da = (nda.min()-1, nda.max()+1)
      ds = da

    else :
      ave, rms = nda.mean(), nda.std()
      da = (ave-3*rms, ave+3*rms)
      ds = (ave-3*rms, ave+6*rms)

    prefix = 'fig-v0-cspad-RadialBkgd'

    gg.plotImageLarge(img, amp_range=da, figsize=(14,12))
    gg.save('%s-%02d-img.png' % (prefix, ntest))

    gg.hist1d(nda, bins=None, amp_range=ds, weights=None, color=None, show_stat=True, log=False, \
           figsize=(6,5), axwin=(0.18, 0.12, 0.78, 0.80), \
           title=None, xlabel='Pixel value', ylabel='Number of pixels', titwin='Pixel value sprctrum')
    gg.save('%s-%02d-his.png' % (prefix, ntest))

    gg.show()

#------------------------------

if __name__ == '__main__' :
    import sys
    ntest = int(sys.argv[1]) if len(sys.argv)>1 else 1
    print 'Test # %d' % ntest
    test(ntest)
    sys.exit('End of test')
 
#------------------------------
#------------------------------
#------------------------------
#------------------------------
