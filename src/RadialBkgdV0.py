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
    orb   = rb.obj_radbins() # returns HBins object for radial bins
    opb   = rb.obj_phibins() # returns HBins object for angular bins
    rad   = rb.pixel_rad()
    irad  = rb.pixel_irad()
    phi0  = rb.pixel_phi0()
    phi   = rb.pixel_phi()
    iphi  = rb.pixel_iphi()
    iseq  = rb.pixel_iseq()
    npix  = rb.npixels_per_bin()
    int   = rb.intensity_per_bin(nda)
    arr   = rb.average_per_bin(nda)
    arr   = rb.average_rad_phi(nda, do_transp=True)
    bkgd  = rb.bkgd_nda(nda)
    bkgd  = rb.bkgd_nda_interpol(nda, method='linear') # method='nearest' 'cubic'
    cdata = rb.subtract_bkgd(nda)
    cdata = rb.subtract_bkgd_interpol(nda, method='linear')


    # Print attributes and n-d arrays
    # -------------------------------
    rb.print_attrs()
    rb.print_ndarrs()

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
from __future__ import print_function
from __future__ import division
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

class RadialBkgd(object) :
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
        self.shapeflat = (self.rad.size,)
        self.rad.shape = self.shapeflat
        self.phi0.shape = self.shapeflat
        self.mask = mask

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

        self.griddata = None
        self.print_ndarr = None


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
        phi1, phi2 = self.pb.limits()
        self.is360 = math.fabs(math.fabs(phi2-phi1)-360) < 1e-3


    def print_attrs(self) :
        print('%s attrbutes:' % self.__class__.__name__)
        print(self.pb.strrange(fmt='Phi bins:  min:%8.1f  max:%8.1f  nbins:%5d'))
        print(self.rb.strrange(fmt='Rad bins:  min:%8.1f  max:%8.1f  nbins:%5d'))


    def print_ndarrs(self) :
        print('%s n-d arrays:' % self.__class__.__name__)
        if self.print_ndarr is None :
            from Detector.GlobalUtils import print_ndarr
            self.print_ndarr = print_ndarr
        self.print_ndarr(self.rad, '  rad')
        self.print_ndarr(self.phi, '  phi')
        self.print_ndarr(self.mask,'  mask')
        #print 'Phi limits: ', phiedges[0], phiedges[-1]


    def obj_radbins(self) :
        """Returns HBins object for radial bins."""
        return self.rb


    def obj_phibins(self) :
        """Returns HBins object for angular bins."""
        return self.pb


    def pixel_rad(self) :
        """Returns 1-d numpy array of pixel radial parameters."""
        return self.rad


    def pixel_irad(self) :
        """Returns 1-d numpy array of pixel radial indexes."""
        return self.irad


    def pixel_phi0(self) :
        """Returns 1-d numpy array of pixel angules in the range [-180,180] degree."""
        return self.phi0


    def pixel_phi(self) :
        """Returns 1-d numpy array of pixel angules in the range [phi_min, phi_min+360] degree."""
        return self.phi


    def pixel_iphi(self) :
        """Returns 1-d numpy array of pixel angular indexes."""
        return self.iphi


    def pixel_iseq(self) :
        """Returns 1-d numpy array of sequentially (in rad and phi) numerated pixel indexes."""
        return self.iseq


    def npixels_per_bin(self) :
        """Returns 1-d numpy array of number of accounted pixels per bin."""
        return self.npix_per_bin


    def _flatten_(self, nda) :
        if len(nda.shape)>1 :
            #nda = nda.flatten()
            nda.shape = self.shapeflat
        return nda


    def intensity_per_bin(self, nda) :
        """Returns 1-d numpy array of total pixel intensity per bin for input array nda."""
        return np.bincount(self.iseq, weights=self._flatten_(nda), minlength=None)


    def average_per_bin(self, nda) :
        """Returns 1-d numpy array of averaged in bin intensity for input array nda."""
        num = self.intensity_per_bin(self._flatten_(nda))
        den = self.npixels_per_bin()
        return divide_protected(num, den, vsub_zero=0)


    def average_rad_phi(self, nda, do_transp=True) :
        """Returns 2-d (rad,phi) numpy array of averaged in bin intensity for input array nda."""
        arr_rphi = self.average_per_bin(self._flatten_(nda))[:-1]
        arr_rphi.shape = (self.pb.nbins(), self.rb.nbins())
        return np.transpose(arr_rphi) if do_transp else arr_rphi


    def bkgd_nda(self, nda) :
        """Returns 1-d numpy array of per-pixel background for input array nda."""
        bin_bkgd = self.average_per_bin(self._flatten_(nda))
        return np.array([bin_bkgd[i] for i in self.iseq])


    def bkgd_nda_interpol(self, nda, method='linear', verb=False) : # 'nearest' 'cubic'
        """Returns 1-d numpy array of per-pixel interpolated background for averaged input data."""

        #if not is360 : raise ValueError('Interpolation works for 360 degree coverage ONLY') 

        if self.griddata is None :
            from scipy.interpolate import griddata
            self.griddata = griddata

        # 1) get values in bin centers
        binv = self.average_rad_phi(self._flatten_(nda), do_transp=False)

        # 2) add values in bin edges
        
        if verb : print('binv.shape: ', binv.shape)
        vrad_a1,  vrad_a2 = binv[0,:], binv[-1,:]
        if self.is360 :
            vrad_a1 = vrad_a2 = 0.5*(binv[0,:] + binv[-1,:]) # [iphi, irad]
        nodea = np.vstack((vrad_a1, binv, vrad_a2))
        
        vang_rmin, vang_rmax = nodea[:,0], nodea[:,-1]
        vang_rmin.shape = vang_rmax.shape = (vang_rmin.size, 1) # it should be 2d for hstack
        val_nodes = np.hstack((vang_rmin, nodea, vang_rmax))
        if verb : print('nodear.shape: ', val_nodes.shape)

        # 3) extend bin-centers by limits        
        bcentsr = self.rb.bincenters()
        bcentsp = self.pb.bincenters()
        blimsr  = self.rb.limits()
        blimsp  = self.pb.limits()

        rad_nodes = np.concatenate(((blimsr[0],), bcentsr, (blimsr[1],)))
        phi_nodes = np.concatenate(((blimsp[0],), bcentsp, (blimsp[1],)))
        if verb : print('rad_nodes.shape', rad_nodes.shape)
        if verb : print('phi_nodes.shape', phi_nodes.shape)

        # 4) make point coordinate and value arrays
        points_rad, points_phi = np.meshgrid(rad_nodes, phi_nodes)
        if verb : print('points_phi.shape', points_phi.shape)
        if verb : print('points_rad.shape', points_rad.shape)
        points = np.array(list(zip(points_phi.flatten(), points_rad.flatten()))) 
        if verb : print('points.shape', points.shape)

        values = val_nodes.flatten()
        if verb : print('values.shape', values.shape)

        # 4) return interpolated data on (phi, rad) grid
        grid_vals = self.griddata(points, values, (self.phi, self.rad), method=method)
        return np.select((self.iseq<self.pb.nbins()*self.rb.nbins(),), (grid_vals,), default=0)


    def subtract_bkgd(self, ndarr) :
        """Returns 1-d numpy array of per-pixel background subtracted input data."""
        nda = self._flatten_(ndarr)
        return nda - self.bkgd_nda(nda)


    def subtract_bkgd_interpol(self, ndarr, method='linear', verb=False) :
        """Returns 1-d numpy array of per-pixel interpolated-background subtracted input data."""
        nda = self._flatten_(ndarr)
        return nda - self.bkgd_nda_interpol(nda, method, verb)

#------------------------------
#------------------------------
#----------- TEST -------------
#------------------------------
#------------------------------

def data_geo(ntest) :
    """Returns test data numpy array and geometry object
    """
    from time import time
    from PSCalib.NDArrIO import save_txt, load_txt
    from PSCalib.GeometryAccess import GeometryAccess

    dir       = '/reg/g/psdm/detector/alignment/cspad/calib-cxi-camera2-2016-02-05'
    #fname_nda = '%s/nda-water-ring-cxij4716-r0022-e000001-CxiDs2-0-Cspad-0-ave.txt' % dir
    fname_nda = '%s/nda-water-ring-cxij4716-r0022-e014636-CxiDs2-0-Cspad-0-ave.txt' % dir
    fname_geo = '%s/calib/CsPad::CalibV1/CxiDs2.0:Cspad.0/geometry/geo-cxi01516-2016-02-18-Ag-behenate-tuned.data' % dir
    #fname_geo = '%s/geo-cxi02416-r0010-2016-03-11.txt' % dir
    fname_gain = '%s/calib/CsPad::CalibV1/CxiDs2.0:Cspad.0/pixel_gain/cxi01516-r0016-2016-02-18-FeKalpha.data' % dir

    # load n-d array with averaged water ring
    arr = load_txt(fname_nda)
    #arr *= load_txt(fname_gain)
    #print_ndarr(arr,'water ring')
    arr.shape = (arr.size,) # (32*185*388,)

    # retrieve geometry
    t0_sec = time()
    geo = GeometryAccess(fname_geo)
    geo.move_geo('CSPAD:V1', 0, 1600, 0, 0)
    geo.move_geo('QUAD:V1', 2, -100, 0, 0)
    #geo.get_geo('QUAD:V1', 3).print_geo()
    print('Time to load geometry %.3f sec from file\n%s' % (time()-t0_sec, fname_geo))

    return arr, geo

#------------------------------

def test01(ntest) :
    """Test for radial 1-d binning of entire image.
    """
    from time import time
    import pyimgalgos.GlobalGraphics as gg
    from PSCalib.GeometryAccess import img_from_pixel_arrays

    arr, geo = data_geo(ntest)

    t0_sec = time()
    iX, iY = geo.get_pixel_coord_indexes()
    X, Y, Z = geo.get_pixel_coords()
    mask = geo.get_pixel_mask(mbits=0o377).flatten() 
    print('Time to retrieve geometry %.3f sec' % (time()-t0_sec))

    t0_sec = time()
    rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=1) # v1
    print('RadialBkgd initialization time %.3f sec' % (time()-t0_sec))

    t0_sec = time()
    nda, title = arr, None
    if   ntest == 1 : nda, title = arr,                   'averaged data'
    elif ntest == 2 : nda, title = rb.pixel_rad(),        'pixel radius value'
    elif ntest == 3 : nda, title = rb.pixel_phi(),        'pixel phi value'
    elif ntest == 4 : nda, title = rb.pixel_irad() + 2,   'pixel radial bin index' 
    elif ntest == 5 : nda, title = rb.pixel_iphi() + 2,   'pixel phi bin index'
    elif ntest == 6 : nda, title = rb.pixel_iseq() + 2,   'pixel sequential (rad and phi) bin index'
    elif ntest == 7 : nda, title = mask,                  'mask'
    elif ntest == 8 : nda, title = rb.bkgd_nda(nda),      'averaged radial background'
    elif ntest == 9 : nda, title = rb.subtract_bkgd(nda) * mask, 'background-subtracted data'

    else :
        t1_sec = time()
        pf = polarization_factor(rb.pixel_rad(), rb.pixel_phi(), 94e3) # Z=94mm
        print('Time to evaluate polarization correction factor %.3f sec' % (time()-t1_sec))

        if   ntest ==10 : nda, title = pf,                    'polarization factor'
        elif ntest ==11 : nda, title = arr * pf,              'polarization-corrected averaged data'
        elif ntest ==12 : nda, title = rb.subtract_bkgd(arr * pf) * mask , 'polarization-corrected background subtracted data'
        elif ntest ==13 : nda, title = rb.bkgd_nda(arr * pf), 'polarization-corrected averaged radial background'
        elif ntest ==14 : nda, title = rb.bkgd_nda_interpol(arr * pf) * mask , 'polarization-corrected interpolated radial background'
        elif ntest ==15 : nda, title = rb.subtract_bkgd_interpol(arr * pf) * mask , 'polarization-corrected interpolated radial background-subtracted data'


        else :
            print('Test %d is not implemented' % ntest) 
            return
        
    print('Get %s n-d array time %.3f sec' % (title, time()-t0_sec))

    img = img_from_pixel_arrays(iX, iY, nda) if not ntest in (21,) else nda[100:300,:]

    da, ds = None, None
    colmap = 'jet' # 'cubehelix' 'cool' 'summer' 'jet' 'winter'
    if ntest in (2,3,4,5,6,7) :
        da = ds = (nda.min()-1., nda.max()+1.)

    if ntest in (12,15) :
        ds = da = (-20, 20)
        colmap = 'gray'

    else :
        ave, rms = nda.mean(), nda.std()
        da = ds = (ave-2*rms, ave+3*rms)

    prefix = 'fig-v01-cspad-RadialBkgd'

    gg.plotImageLarge(img, amp_range=da, figsize=(14,12), title=title, cmap=colmap)
    gg.save('%s-%02d-img.png' % (prefix, ntest))

    gg.hist1d(nda, bins=None, amp_range=ds, weights=None, color=None, show_stat=True, log=False, \
           figsize=(6,5), axwin=(0.18, 0.12, 0.78, 0.80), \
           title=None, xlabel='Pixel value', ylabel='Number of pixels', titwin=title)
    gg.save('%s-%02d-his.png' % (prefix, ntest))

    gg.show()

    print('End of test for %s' % title)    

#------------------------------

def test02(ntest) :
    """Test for 2-d (default) binning of the rad-phi range of entire image
    """
    #from Detector.GlobalUtils import print_ndarr
    from time import time
    import pyimgalgos.GlobalGraphics as gg
    from PSCalib.GeometryAccess import img_from_pixel_arrays

    arr, geo = data_geo(ntest)

    iX, iY = geo.get_pixel_coord_indexes()
    X, Y, Z = geo.get_pixel_coords()
    mask = geo.get_pixel_mask(mbits=0o377).flatten() 

    t0_sec = time()
    rb = RadialBkgd(X, Y, mask) # v0
    #rb = RadialBkgd(X, Y, mask, nradbins=500) # , nphibins=8, phiedges=(-20, 240), radedges=(10000,80000))
    print('RadialBkgd initialization time %.3f sec' % (time()-t0_sec))

    #print 'npixels_per_bin:',   rb.npixels_per_bin()
    #print 'intensity_per_bin:', rb.intensity_per_bin(arr)
    #print 'average_per_bin:',   rb.average_per_bin(arr)

    t0_sec = time()
    nda, title = arr, None
    if   ntest == 21 : nda, title = arr,                   'averaged data'
    elif ntest == 22 : nda, title = rb.pixel_rad(),        'pixel radius value'
    elif ntest == 23 : nda, title = rb.pixel_phi(),        'pixel phi value'
    elif ntest == 24 : nda, title = rb.pixel_irad() + 2,   'pixel radial bin index' 
    elif ntest == 25 : nda, title = rb.pixel_iphi() + 2,   'pixel phi bin index'
    elif ntest == 26 : nda, title = rb.pixel_iseq() + 2,   'pixel sequential (rad and phi) bin index'
    elif ntest == 27 : nda, title = mask,                  'mask'
    elif ntest == 28 : nda, title = rb.bkgd_nda(nda),      'averaged radial background'
    elif ntest == 29 : nda, title = rb.subtract_bkgd(nda) * mask, 'background-subtracted data'
    elif ntest == 30 : nda, title = rb.average_rad_phi(nda),'r-phi'
    elif ntest == 31 : nda, title = rb.bkgd_nda_interpol(nda), 'averaged radial interpolated background'
    elif ntest == 32 : nda, title = rb.subtract_bkgd_interpol(nda, method='linear', verb=True) * mask, 'interpol-background-subtracted data'
    else :
        print('Test %d is not implemented' % ntest) 
        return

    print('Get %s n-d array time %.3f sec' % (title, time()-t0_sec))

    img = img_from_pixel_arrays(iX, iY, nda) if not ntest in (30,) else nda # [100:300,:]

    colmap = 'jet' # 'cubehelix' 'cool' 'summer' 'jet' 'winter' 'gray'

    da = (nda.min()-1, nda.max()+1)
    ds = da

    if ntest in (21,28,29,30,31) :
        ave, rms = nda.mean(), nda.std()
        da = ds = (ave-2*rms, ave+3*rms)

    elif ntest in (32,) : 
        colmap = 'gray'
        ds = da = (-20, 20)

    prefix = 'fig-v02-cspad-RadialBkgd'

    gg.plotImageLarge(img, amp_range=da, figsize=(14,12), title=title, cmap=colmap)
    gg.save('%s-%02d-img.png' % (prefix, ntest))

    gg.hist1d(nda, bins=None, amp_range=ds, weights=None, color=None, show_stat=True, log=False, \
           figsize=(6,5), axwin=(0.18, 0.12, 0.78, 0.80), \
           title=None, xlabel='Pixel value', ylabel='Number of pixels', titwin=title)
    gg.save('%s-%02d-his.png' % (prefix, ntest))

    gg.show()

    print('End of test for %s' % title)    

#------------------------------

def test03(ntest) :
    """Test for 2-d binning of the restricted rad-phi range of entire image
    """
    from time import time
    import pyimgalgos.GlobalGraphics as gg
    from PSCalib.GeometryAccess import img_from_pixel_arrays

    arr, geo = data_geo(ntest)

    iX, iY = geo.get_pixel_coord_indexes()
    X, Y, Z = geo.get_pixel_coords()
    mask = geo.get_pixel_mask(mbits=0o377).flatten() 

    t0_sec = time()

    rb = RadialBkgd(X, Y, mask, nradbins=200, nphibins=32, phiedges=(-20, 240), radedges=(10000,80000)) if ntest in (51,52)\
    else RadialBkgd(X, Y, mask, nradbins=  5, nphibins= 8, phiedges=(-20, 240), radedges=(10000,80000))
    #rb = RadialBkgd(X, Y, mask, nradbins=3, nphibins=8, phiedges=(240, -20), radedges=(80000,10000)) # v3

    print('RadialBkgd initialization time %.3f sec' % (time()-t0_sec))

    #print 'npixels_per_bin:',   rb.npixels_per_bin()
    #print 'intensity_per_bin:', rb.intensity_per_bin(arr)
    #print 'average_per_bin:',   rb.average_per_bin(arr)

    t0_sec = time()
    nda, title = arr, None
    if   ntest == 41 : nda, title = arr,                   'averaged data'
    elif ntest == 42 : nda, title = rb.pixel_rad(),        'pixel radius value'
    elif ntest == 43 : nda, title = rb.pixel_phi(),        'pixel phi value'
    elif ntest == 44 : nda, title = rb.pixel_irad() + 2,   'pixel radial bin index' 
    elif ntest == 45 : nda, title = rb.pixel_iphi() + 2,   'pixel phi bin index'
    elif ntest == 46 : nda, title = rb.pixel_iseq() + 2,   'pixel sequential (rad and phi) bin index'
    elif ntest == 47 : nda, title = mask,                  'mask'
    elif ntest == 48 : nda, title = rb.bkgd_nda(nda),      'averaged radial background'
    elif ntest == 49 : nda, title = rb.subtract_bkgd(nda) * mask, 'background-subtracted data'
    elif ntest == 50 : nda, title = rb.average_rad_phi(nda),'r-phi'
    elif ntest == 51 : nda, title = rb.bkgd_nda_interpol(nda), 'averaged radial interpolated background'
    elif ntest == 52 : nda, title = rb.subtract_bkgd_interpol(nda) * mask, 'interpol-background-subtracted data'
    else :
        print('Test %d is not implemented' % ntest) 
        return

    print('Get %s n-d array time %.3f sec' % (title, time()-t0_sec))

    img = img_from_pixel_arrays(iX, iY, nda) if not ntest in (50,) else nda # [100:300,:]

    colmap = 'jet' # 'cubehelix' 'cool' 'summer' 'jet' 'winter' 'gray'

    da = (nda.min()-1, nda.max()+1)
    ds = da

    if ntest in (41,48,49,50,51) :
        ave, rms = nda.mean(), nda.std()
        da = ds = (ave-2*rms, ave+3*rms)

    elif ntest in (52,) : 
        colmap = 'gray'
        ds = da = (-20, 20)

    prefix = 'fig-v03-cspad-RadialBkgd'

    gg.plotImageLarge(img, amp_range=da, figsize=(14,12), title=title, cmap=colmap)
    gg.save('%s-%02d-img.png' % (prefix, ntest))

    gg.hist1d(nda, bins=None, amp_range=ds, weights=None, color=None, show_stat=True, log=False, \
           figsize=(6,5), axwin=(0.18, 0.12, 0.78, 0.80), \
           title=None, xlabel='Pixel value', ylabel='Number of pixels', titwin=title)
    gg.save('%s-%02d-his.png' % (prefix, ntest))

    gg.show()

    print('End of test for %s' % title)    

#------------------------------

if __name__ == '__main__' :
    import sys
    ntest = int(sys.argv[1]) if len(sys.argv)>1 else 1
    print('Test # %d' % ntest)
    if   ntest<20 : test01(ntest)
    elif ntest<40 : test02(ntest)
    elif ntest<60 : test03(ntest)
    else : print('Test %d is not implemented' % ntest)     
    #sys.exit('End of test')
 
#------------------------------
#------------------------------
