#!/usr/bin/env python

from __future__ import print_function
from pyimgalgos.RadialBkgd import RadialBkgd, polarization_factor

import math
import numpy as np
#from pyimgalgos.HBins import HBins

from time import time
#from Detector.GlobalUtils import print_ndarr
import pyimgalgos.GlobalGraphics as gg
from PSCalib.GeometryAccess import GeometryAccess, img_from_pixel_arrays
from PSCalib.NDArrIO import save_txt, load_txt

#from pyimgalgos.NDArrGenerators import random_standard

#------------------------------

def test(ntest) :

    #arr = random_standard(shape=(40,60), mu=200, sigma=25)

    #import psana
    #ds  = psana.DataSource('exp=cxij4716:run=22')
    #det = psana.Detector('CxiDs2.0:Cspad.0', ds.env())

    prefix = 'fig-v00-cspad-RadialBkgd'

    # data v1
    #dir       = '/reg/g/psdm/detector/alignment/cspad/calib-cxi-camera2-2016-02-05'
    ##fname_nda = '%s/nda-water-ring-cxij4716-r0022-e000001-CxiDs2-0-Cspad-0-ave.txt' % dir
    ##fname_nda = '%s/nda-water-ring-cxij4716-r0022-e001000-CxiDs2-0-Cspad-0-ave.txt' % dir
    #fname_nda = '%s/nda-water-ring-cxij4716-r0022-e014636-CxiDs2-0-Cspad-0-ave.txt' % dir
    #fname_geo = '%s/geo-cxi02416-r0010-2016-03-11.txt' % dir


    # data v2
    #dir = '/reg/g/psdm/detector/alignment/cspad/calib-cxi-camera2-2015-01-20/2016-03-21-cxih9615-Oleksandr'
    ##fname_nda = '%s/nda-t4-cxih9615-r0039-e030606-CxiDs2-0-Cspad-0-ave.txt' % dir
    #fname_nda = '%s/nda-t4-cxih9615-r0051-e117611-CxiDs2-0-Cspad-0-ave.txt' % dir
    #fname_geo = '%s/geo-cxih9615-r48-51-65-v1.data' % dir


    # data v3
    dir = '/reg/g/psdm/detector/alignment/cspad/calib-cxi-camera2-2015-01-20/'
    ##fname_geo = '%s/geo-cxii5615-r24-42-v1.data' % dir # initial geometry
    fname_geo = '%s/geo-cxij4915-r25-v1.data' % dir # aligned geometry using water ring
    fname_nda = '%s/cspad-ndarr-ave-cxij4915-r0025-e106262-water-ring.txt' % dir
    ##fname_nda = '%s/cspad-ndarr-ave-cxij4915-r0026-e094241-water-ring.txt' % dir
    prefix = 'fig-v03-cspad-RadialBkgd'


    # data v4
    #dir = '/reg/g/psdm/detector/alignment/cspad/calib-cxi-camera2-2015-01-20/'
    #fname_geo = '%s/geo-cxih9615-r48-51-65-v1.data' % dir
    #fname_nda = '%s/cspad-ndarr-max-cxih9615-r48-51-65-protein-rings.txt' % dir
    #prefix = 'fig-v04-cspad-RadialBkgd'


    # data v5
    #dir = '/reg/g/psdm/detector/alignment/cspad/calib-cxi-camera2-2015-01-20/2015-11-06-cxif5315-geo-tuning'
    #fname_geo = '%s/geo-cxif5315-r0169-2015-11-06-tuned.data' % dir # aligned geometry using water ring
    #fname_nda = '%s/nda-ave-cxif5315-r0169-85188ev-data.txt' % dir
    #prefix = 'fig-v05-nbin32-cspad-RadialBkgd'


    # data v6
    #dir = '/reg/g/psdm/detector/alignment/cspad/calib-cxi-camera2-2015-01-20/2015-11-06-cxif5315-geo-tuning'
    #fname_geo = '%s/geo-cxif5315-r0169-2015-11-06-tuned.data' % dir # aligned geometry using water ring
    #fname_nda = 'nda-cxi00516-r0014-e000100-CxiDs2-0-Cspad-0-ave.txt'
    ##fname_nda = '%s/nda-max-cxif5315-r0169-85188ev-data.txt' % dir
    #prefix = 'fig-v06-cspad-RadialBkgd'


    # load n-d array with averaged water ring
    arr = load_txt(fname_nda)
    arr.shape = (arr.size,) # (32*185*388,)

    # retrieve geometry
    t0_sec = time()
    geo = GeometryAccess(fname_geo)

    #geo.move_geo('QUAD:V1', 0, 0, 200, 0)
    #geo.move_geo('QUAD:V1', 1, 0, 200, 0)

    #geo.move_geo('CSPAD:V1', 0, 1600, 0, 0)
    #geo.move_geo('QUAD:V1', 2, -100, 0, 0)
    #geo.move_geo('QUAD:V1', 3, -300, -100, 0)
    #geo.move_geo('QUAD:V1', 3, 500, 500, 0)
    #geo.tilt_geo('QUAD:V1', 3, 0, 0, 1)
    #geo.move_geo('QUAD:V1', 3, -500, -200, 0)
    #geo.get_geo('QUAD:V1', 3).print_geo()

    # data v3
    #geo.move_geo('CSPAD:V1', 0, -2200, -600, 0) # for geo-cxii5615-r24-42-v1.data

    iX, iY = geo.get_pixel_coord_indexes()
    X, Y, Z = geo.get_pixel_coords()
    mask = geo.get_pixel_mask(mbits=0o377).flatten() 

    print('Time to retrieve geometry %.3f sec' % (time()-t0_sec))

    t0_sec = time()
    #rb = RadialBkgd(X, Y, mask) # v0
    #rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=32) # v1
    #rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=8) # v1
    rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=1) # v1
    #rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=200) # v5
    #rb = RadialBkgd(X, Y, mask, nradbins=500, nphibins=8, phiedges=(-20, 240), radedges=(10000,80000)) # v2
    #rb = RadialBkgd(X, Y, mask, nradbins=3, nphibins=8, phiedges=(240, -20), radedges=(80000,10000)) # v3
    #rb = RadialBkgd(X, Y, mask, nradbins=3, nphibins=8, phiedges=(-20, 240), radedges=(10000,80000))

    print('RadialBkgd initialization time %.3f sec' % (time()-t0_sec))

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
    elif ntest == 21: nda, title = rb.average_rad_phi(nda),'r-phi'

    else :
        t1_sec = time()
        #pf = polarization_factor(rb.pixel_rad(), rb.pixel_phi(), 94e3) # Z=94mm
        pf = polarization_factor(rb.pixel_rad(), rb.pixel_phi(), 91.33e3) # Z=913.3mm
        print('Time to evaluate polarization correction factor %.3f sec' % (time()-t1_sec))

        if   ntest ==10 : nda, title = pf,                    'polarization factor'
        elif ntest ==11 : nda, title = arr * pf,              'polarization-corrected averaged data'
        elif ntest ==12 : nda, title = rb.subtract_bkgd(arr * pf) * mask , 'polarization-corrected radial background-subtracted data'
        elif ntest ==13 : nda, title = rb.bkgd_nda(arr * pf), 'polarization-corrected background'
        elif ntest ==14 : nda, title = rb.bkgd_nda_interpol(arr * pf) * mask , 'polarization-corrected interpolated radial background'
        elif ntest ==15 : nda, title = rb.subtract_bkgd_interpol(arr * pf) * mask , 'polarization-corrected interpolated radial background-subtracted data'

    print('Get %s n-d array time %.3f sec' % (title, time()-t0_sec))

    img = img_from_pixel_arrays(iX, iY, nda) if not ntest in (21,) else nda[100:300,:]

    da, ds = None, None
    colmap = 'jet' # 'cubehelix' 'cool' 'summer' 'jet' 'winter'
    if ntest in (2,3,4,5,6,7) :
        ds = da = (nda.min()-1, nda.max()+1)

    elif ntest in (12,15) :
        #ds = da = (-1600, 2400)
        ds = da = (-20, 20)
        colmap = 'gray'

    else :
        ave, rms = nda.mean(), nda.std()
        da = ds = (ave-2*rms, ave+3*rms)

    save_txt('nda-%s-%02d.txt' % (prefix, ntest), nda, fmt='%.2f', verbos=True)

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
    test(ntest)
    #sys.exit('End of test')
 
#------------------------------
#------------------------------
#------------------------------
#------------------------------
