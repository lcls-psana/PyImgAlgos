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
    phi0         = rb.pixel_phi0()
    phi          = rb.pixel_phi()
    iphi         = rb.pixel_iphi()
    iseq         = rb.pixel_iseq()
    npix_per_bin = rb.npixels_per_bin()
    int_per_bin  = rb.intensity_per_bin(nda)
    bkgd         = rb.bkgd_nda(nda)
    res          = rb.subtract_bkgd(nda)

    # Global methods
    # --------------

    from pyimgalgos.RadialBkgd import polarization_factor

    polf = polarization_factor(rad, phi, z)
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


def polarization_factor(rad, phi_deg, z) :
    """Returns per-pixel polarization factors, assuming that detector is perpendicular to Z.
    """
    phi = np.deg2rad(phi_deg)
    ones = np.ones_like(rad)
    theta = np.arctan2(rad, z)
    sxc = np.sin(theta)*np.cos(phi)
    pol = 1 - sxc*sxc
    return divide_protected(ones, pol, vsub_zero=0)

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
        self.rad, self.phi0 = cart2polar(xarr, yarr)
        shape = (self.rad.size,)
        self.rad.shape = shape
        self.phi0.shape = shape

        phimin = min(phiedges[0], phiedges[-1])

        self.phi = np.select((self.phi0<phimin, self.phi0>=phimin), (self.phi0+360.,self.phi0))

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

        if mask is not None : cond *= mask.flatten()

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
        #print 'Phi limits: ', phiedges[0], phiedges[-1]


    def pixel_obj_radbins(self) :
        """Returns HBins object for radial bins."""
        return self.rb


    def pixel_obj_phibins(self) :
        """Returns HBins object for angular bins."""
        return self.pb


    def pixel_rad(self) :
        """Returns numpy array of pixel radial parameters."""
        return self.rad


    def pixel_irad(self) :
        """Returns numpy array of pixel radial indexes."""
        return self.irad


    def pixel_phi0(self) :
        """Returns numpy array of pixel angules in the range [-180,180] degree."""
        return self.phi0


    def pixel_phi(self) :
        """Returns numpy array of pixel angules in the range [phi_min, phi_min+360] degree."""
        return self.phi


    def pixel_iphi(self) :
        """Returns numpy array of pixel angular indexes."""
        return self.iphi


    def pixel_iseq(self) :
        """Returns numpy array of sequentially (in rad and phi) numerated pixel indexes."""
        return self.iseq


    def npixels_per_bin(self) :
        """Returns numpy array of number of accounted pixels per bin."""
        return self.npix_per_bin


    def intensity_per_bin(self, nda) :
        """Returns numpy array of total pixel intensity per bin for input array nda."""
        return np.bincount(self.iseq, weights=nda.flatten(), minlength=None)


    def average_per_bin(self, nda) :
        """Returns numpy array of averaged in bin intensity for input array nda."""
        num = self.intensity_per_bin(nda)
        den = self.npixels_per_bin()
        return divide_protected(num, den, vsub_zero=0)


    def bkgd_nda(self, nda) :
        """Returns numpy array of per-pixel background for input array nda."""
        bin_bkgd = self.average_per_bin(nda)
        return np.array([bin_bkgd[i] for i in self.iseq])


    def subtract_bkgd(self, nda) :
        """Returns numpy array of background subtracted input array nda."""
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
    #fname_nda = '%s/nda-water-ring-cxij4716-r0022-e000001-CxiDs2-0-Cspad-0-ave.txt' % dir
    #fname_nda = '%s/nda-water-ring-cxij4716-r0022-e001000-CxiDs2-0-Cspad-0-ave.txt' % dir
    fname_nda = '%s/nda-water-ring-cxij4716-r0022-e014636-CxiDs2-0-Cspad-0-ave.txt' % dir
    fname_geo = '%s/calib/CsPad::CalibV1/CxiDs2.0:Cspad.0/geometry/geo-cxi01516-2016-02-18-Ag-behenate-tuned.data' % dir

    # load n-d array with averaged water ring
    arr = load_txt(fname_nda)
    #print_ndarr(arr,'water ring')
    arr.shape = (arr.size,) # (32*185*388,)

    # retrieve geometry
    t0_sec = time()
    geo = GeometryAccess(fname_geo)
    iX, iY = geo.get_pixel_coord_indexes()
    X, Y, Z = geo.get_pixel_coords()
    mask = geo.get_pixel_mask(mbits=0377).flatten() 
    print 'Time to retrieve geometry %.3f sec' % (time()-t0_sec)

    t0_sec = time()
    #rb = RadialBkgd(X, Y, mask) # v0
    rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=8) # v1
    #rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=1) # v1
    #rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=200) # v5
    #rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=8, phiedges=(-20, 240), radedges=(10000,80000)) # v2
    #rb = RadialBkgd(X, Y, mask, nradbins=3, nphibins=8, phiedges=(240, -20), radedges=(80000,10000)) # v3
    #rb = RadialBkgd(X, Y, mask, nradbins=3, nphibins=8, phiedges=(-20, 240), radedges=(10000,80000))

    print 'RadialBkgd initialization time %.3f sec' % (time()-t0_sec)

    #print 'npixels_per_bin:',   rb.npixels_per_bin()
    #print 'intensity_per_bin:', rb.intensity_per_bin(arr)
    #print 'average_per_bin:',   rb.average_per_bin(arr)

    t0_sec = time()
    nda, title = arr, 'averaged data'
    if   ntest == 2 : nda, title = rb.pixel_rad(),        'pixel radius value'
    elif ntest == 3 : nda, title = rb.pixel_phi(),        'pixel phi value'
    elif ntest == 4 : nda, title = rb.pixel_irad() + 2,   'pixel radial bin index' 
    elif ntest == 5 : nda, title = rb.pixel_iphi() + 2,   'pixel phi bin index'
    elif ntest == 6 : nda, title = rb.pixel_iseq() + 2,   'pixel sequential (inr and phi) bin index'
    elif ntest == 7 : nda, title = mask,                  'mask'
    elif ntest == 8 : nda, title = rb.bkgd_nda(nda),      'averaged radial background'
    elif ntest == 9 : nda, title = rb.subtract_bkgd(nda) * mask, 'background-subtracted data'
    else :
        t1_sec = time()
        pf = polarization_factor(rb.pixel_rad(), rb.pixel_phi(), 94e3) # Z=94mm
        print 'Time to evaluate polarization correction factor %.3f sec' % (time()-t1_sec)

        if   ntest ==10 : nda, title = pf,                    'polarization factor'
        elif ntest ==11 : nda, title = arr * pf,              'polarization-corrected averaged data'
        elif ntest ==12 : nda, title = rb.subtract_bkgd(arr * pf) * mask , 'polarization-corrected radial background'
        elif ntest ==13 : nda, title = rb.bkgd_nda(arr * pf), 'polarization-corrected background-subtracted data'

    print 'Get %s n-d array time %.3f sec' % (title, time()-t0_sec)

    img = img_from_pixel_arrays(iX, iY, nda)

    da, ds = None, None
    colmap = 'jet' # 'cubehelix' 'cool' 'summer' 'jet' 'winter'
    if ntest in (2,3,4,5,6,7) :
        da = (nda.min()-1, nda.max()+1)
        ds = da

    if ntest in (12,) :
        ds = da = (-20, 20)
        colmap = 'gray'
    else :
        ave, rms = nda.mean(), nda.std()
        da = ds = (ave-2*rms, ave+3*rms)

    prefix = 'fig-v6-cspad-RadialBkgd'

    gg.plotImageLarge(img, amp_range=da, figsize=(14,12), title=title, cmap=colmap)
    gg.save('%s-%02d-img.png' % (prefix, ntest))

    gg.hist1d(nda, bins=None, amp_range=ds, weights=None, color=None, show_stat=True, log=False, \
           figsize=(6,5), axwin=(0.18, 0.12, 0.78, 0.80), \
           title=None, xlabel='Pixel value', ylabel='Number of pixels', titwin=title)
    gg.save('%s-%02d-his.png' % (prefix, ntest))

    gg.show()

    print 'End of test for %s' % title    

#------------------------------

if __name__ == '__main__' :
    import sys
    ntest = int(sys.argv[1]) if len(sys.argv)>1 else 1
    print 'Test # %d' % ntest
    test(ntest)
    #sys.exit('End of test')
 
#------------------------------
#------------------------------
#------------------------------
#------------------------------
