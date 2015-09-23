#!/usr/bin/env python

##-----------------------------
"""A set of methods to evaluate angles in fiber diffraction experiments.

Usage::
    # Import
    from pyimgalgos.FiberAngles import fraser, calc_phi, calc_beta, funcy

    # Fraser transformation:
    s12rot, s3rot, reciparr = fraser(arr, beta_deg, L, center=None, oshape=(1500,1500))

    # Evaluation of fiber tilt angles beta phi (in the image plane) (in transverse to image plane).
    phi  = calc_phi (x1, y1, x2, y2, dist)
    beta = calc_beta(x1, y1, phi, dist)

    Fit functions
    yarr = funcy(xarr, phi_deg, bet_deg)
    yarr2 = funcy2(xarr, a, b, c)

    # Commands to test in the release directory: 
    python ./pyimgalgos/src/FiberAngles.py <test-id>
    # where
    # <test-id> = 1 - test of the Fraser transformation
    # <test-id> = 2 - test of the phi angle
    # <test-id> = 3 - test of the beta angle 
    
This software was developed for the SIT project.
If you use all or part of it, please give an appropriate acknowledgment.

@version $Id$

@author Mikhail S. Dubrovin
"""

##-----------------------------
import numpy as np
import math # pi, sin, sqrt, ceil, floor, fabs
##-----------------------------

class Storage :
    """Storage class for local data exchange between methods.
    """
    def __init__(self) :
        pass

#------------------------------
sp = Storage() # singleton
##-----------------------------

def _fillimg(r,c,a) :
    """This method is used in fraser(...), is called in map(_fillimg, irows, icols, arr) and serves to fill image.
    """
    sp.image[r,c] += a
    sp.count[r,c] += 1
    return r
    
##-----------------------------

def fraser(arr, beta_deg, L, center=None, oshape=(1500,1500)) :
    """Do fraser correction for an array at angle beta and distance L (given in
    units of pixels (110um) - wavelength not considered, angle is in degrees
    example fraser(array,10,909); (10 degrees at 100mm distance)

    ASSUMPTION:
    1. 2-d array center corresponds to image center (0,0) by default
    @param arr      - [in] 2-d image array
    @param beta_deg - [in] angle beta in degrees
    @param L        - [in] distance from sample to detector given in units of pixels (110um)
    @param center   - [in] center corresponds to image center (0,0) by default
    @param oshape   - [in] ouitput image shape
    """

    sizex = arr.shape[0]
    sizey = arr.shape[1]

    xc, yc = center if center is not None else (sizex/2, sizey/2) 

    xarr = np.arange(math.floor(-xc), math.floor(sizex-xc))
    yarr = np.arange(math.floor(-yc), math.floor(sizey-yc))

    x,y = np.meshgrid(yarr, xarr) ### SWAPPED yarr, xarr to keep correct shape for grids

    d = np.sqrt(x*x+y*y+L*L)
    s1 = x/d
    s2 = L/d - 1
    s3 = y/d

    cosbeta = math.cos(math.radians(beta_deg))
    sinbeta = math.sin(math.radians(beta_deg))

    s1rot = s1
    s2rot = s2 * cosbeta - s3 * sinbeta
    s3rot = s2 * sinbeta + s3 * cosbeta

    s12rot = np.sqrt(np.square(s1rot) + np.square(s2rot))
    s12rot[:,1:math.floor(sizex-xc)] *= -1

    s12rot = np.ceil(s12rot * math.floor(sizex-xc))
    s3rot  = np.ceil(s3rot  * math.floor(sizey-yc))

    orows, orows1 = oshape[0], oshape[0] - 1
    ocols, ocols1 = oshape[1], oshape[1] - 1
    
    icols = np.array(s12rot + math.ceil(ocols/2), dtype=np.int)
    irows = np.array(s3rot  + math.ceil(orows/2), dtype=np.int)

    irows = np.select([irows<0, irows>orows1], [0,orows1], default=irows)
    icols = np.select([icols<0, icols>ocols1], [0,ocols1], default=icols)

    #reciparr = np.zeros(oshape, dtype=arr.dtype)
    #counts   = np.zeros(oshape, dtype=np.int)
    #reciparr[irows, icols] = arr

    sp.image = np.zeros(oshape, dtype=arr.dtype)
    sp.count = np.zeros(oshape, dtype=np.int)

    unused_lst = map(_fillimg, irows, icols, arr)

    #print 'arr.shape: ', arr.shape
    #print 's3rot.shape: ', s3rot.shape
    #print 's12rot.shape: ', s12rot.shape
    #print 'reciparr.shape: ', reciparr.shape
    #print 'count min=%d, max=%d' % (sp.count.min(), sp.count.max())    

    countpro = np.select([sp.count<1], [-1], default=sp.count)
    reciparr = np.select([countpro>0], [sp.image/countpro], default=0)

    return s12rot, s3rot, reciparr

##-----------------------------

def calc_phi(x1pix, y1pix, x2pix, y2pix, dist) :
    """Evaluates fiber phi angle [rad] for two peaks in equatorial region
    @param x1pix - [in] x coordinate of the 1st point
    @param y1pix - [in] y coordinate of the 1st point 
    @param x1pix - [in] x coordinate of the 2nd point 
    @param y1pix - [in] y coordinate of the 2nd point 
    @param dist  - [in] distance from sample to detector
    """	
    x1 = x1pix / dist
    y1 = y1pix / dist
    x2 = x2pix / dist
    y2 = y2pix / dist	
    d1 = math.sqrt(x1*x1 + y1*y1 + 1.) - 1.
    d2 = math.sqrt(x2*x2 + y2*y2 + 1.) - 1.
    return math.atan(-(y2*d1 - y1*d2) / (x2*d1 - x1*d2)) 

##-----------------------------

def calc_beta(xpix, ypix, phi, dist) :
    """Evaluates fiber beta angle [rad] for two peaks in equatorial region
    @param xpix - [in] x coordinate of the point
    @param ypix - [in] y coordinate of the point 
    @param phi  - [in] fiber tilt angle [rad] if the detector plane
    @param dist - [in] distance from sample to detector
    """	
    x1 = xpix / dist
    y1 = ypix / dist
    fac = 1. - math.sqrt(x1*x1 + y1*y1 + 1.)
    return -math.atan((y1*math.cos(phi) + x1*math.sin(phi)) / fac)

##-----------------------------

def funcy(x, phi_deg, bet_deg) :
    """Function for parameterization of y(x, phi, beta)
       of peaks in mediane plane for fiber diffraction
       ATTENTION!: curve_fit assumes that x and returned y are numpy arrays.
    """
    phi, bet = math.radians(phi_deg), math.radians(bet_deg)
    sb, cb = math.sin(bet), math.cos(bet)
    t = sb/cb if cb else None
    s, c  = math.sin(phi), math.cos(phi)
    D = c*c - t*t
    if D==0 : return 10
    B = c*(x*s+t)/D
    C = (2*t*x*s + x*x*(s-t)*(s+t))/D
    sqarg = B*B-C
    #print 'sqarg: ', sqarg
    #if sqarg<0 : sqarg=0
    sqapro = np.select([sqarg>0, sqarg<0], [sqarg, 0], default=sqarg)
    sign = 1 if bet>0 else -1
    return -B + sign*np.sqrt(sqapro)

##-----------------------------

def funcy2(x, a, b, c) :
    """Quadratic polynomial function to test curve_fit.
    """    
    return a*x*x + b*x + c

##-----------------------------
##---------- TESTS ------------
##-----------------------------

def test_plot_phi() :
    print """Test plot for phi angle"""

    import pyimgalgos.GlobalGraphics as gg

    xarr = np.linspace(-2,2,50)
    tet = -12
    y0 = [funcy(x,   0, tet) for x in xarr]
    y1 = [funcy(x,  -5, tet) for x in xarr]
    y2 = [funcy(x,  -6, tet) for x in xarr]
    y3 = [funcy(x,  -7, tet) for x in xarr]
    y4 = [funcy(x, -10, tet) for x in xarr]
    y5 = [funcy(x,  10, tet) for x in xarr]
    
    fig1, ax1 = gg.plotGraph(xarr, y0, figsize=(10,5), window=(0.15, 0.10, 0.78, 0.80))
    ax1.plot(xarr, y1,'r-')
    ax1.plot(xarr, y2,'y-')
    ax1.plot(xarr, y3,'k-')
    ax1.plot(xarr, y4,'m-')
    ax1.plot(xarr, y5,'g.')
    ax1.set_xlabel('x', fontsize=14)
    ax1.set_ylabel('y', fontsize=14)
    ax1.set_title('tet=-12, phi=10,0,-5,-6,-7,-10', color='k', fontsize=20)

    #gg.savefig('variation-phi.png')
    gg.show()

##-----------------------------

def test_plot_beta() :
    print """Test plot for beta angle"""

    import pyimgalgos.GlobalGraphics as gg

    xarr = np.linspace(-2,2,50)
    phi = 0
    y0 = [funcy(x, phi,   0) for x in xarr]
    y1 = [funcy(x, phi,  -2) for x in xarr]
    y2 = [funcy(x, phi,  -5) for x in xarr]
    y3 = [funcy(x, phi,  -6) for x in xarr]
    y4 = [funcy(x, phi,  -7) for x in xarr]
    y5 = [funcy(x, phi, -10) for x in xarr]
    y6 = [funcy(x, phi,   2) for x in xarr]
    y7 = [funcy(x, phi,   5) for x in xarr]
    y8 = [funcy(x, phi,  10) for x in xarr]
    
    fig2, ax2 = gg.plotGraph(xarr, y0, figsize=(10,5), window=(0.15, 0.10, 0.78, 0.80)) 
    ax2.plot(xarr, y1,'r-')
    ax2.plot(xarr, y2,'y-')
    ax2.plot(xarr, y3,'b-')
    ax2.plot(xarr, y4,'m-')
    ax2.plot(xarr, y5,'g-')
    ax2.plot(xarr, y6,'r.')
    ax2.plot(xarr, y7,'g.')
    ax2.plot(xarr, y8,'b.')
    ax2.set_xlabel('x', fontsize=14)
    ax2.set_ylabel('y', fontsize=14)
    ax2.set_title('phi=0, theta=10, 5, 2, 0,-2,-5,-6,-7,-10', color='k', fontsize=20)

    #gg.savefig('variation-theta.png')
    gg.show()

##-----------------------------

def test_fraser() :
    print """Test fraser transformation"""

    from pyimgalgos.GlobalGraphics import fig_axes, plot_img, store
    import matplotlib.pyplot as plt

    fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/plots/cspad-cxif5315-0169-000079.npy'

    img2d = np.load(fname)

    s12, s3, recimg = fraser(img2d, 10, 1000, center=None)

    store.fig, store.axim, store.axcb = fig_axes() # if not do_plot else (None, None, None)
    plot_img(recimg, mode='do not hold', amin=0, amax=100)
    #plot_img(img2d, mode='do not hold', amin=0, amax=200)
    plt.ioff() # hold control on show() after the last image
    plt.show()

#------------------------------

if __name__ == "__main__" :

    import sys

    if len(sys.argv)<2 :
        print 'For specific test use command:\n> %s <test-id-string>' % sys.argv[0]
        print 'Default test: test_fraser()'
        test_fraser()
        sys.exit('Default test is completed')

    print 'Test: %s' % sys.argv[1]
    if   sys.argv[1]=='1' : test_fraser()
    elif sys.argv[1]=='2' : test_plot_phi()
    elif sys.argv[1]=='3' : test_plot_beta()
    else :
        print 'Default test: test_fraser()'
        test_fraser()
    sys.exit('Test %s is completed' % sys.argv[1])
    
##-----------------------------
##-----------------------------
##-----------------------------

