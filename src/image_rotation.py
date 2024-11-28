#!/usr/bin/env python

import os
import sys
import math
import numpy as np

class Store(object) :
    """Store of shared parameters.
    """
    def __init__(self) :
        print('In Store')

sp = Store()

#------------------------------

def _fillimg(r,c,a) :
    """This method is used in fraser(...), is called in map(_fillimg, irows, icols, arr) and serves to fill image.
    """
    sp.image[r,c] += a
    sp.count[r,c] += 1
    return r

#------------------------------

def image_rotation(arr, phi_deg, center=None, oshape=(3000,3000)) :
    """Rotate 2-d image arr by phi_deg degree
       Example: image_rotation(arr2d,10);

       ASSUMPTION:
       1. by default 2-d arr image center corresponds to (x,y) origin
       - arr      - [in] 2-d image array
       - phi_deg  - [in] angle phi in degrees
       - center   - [in] center (row,column) location on image, which will be used as (x,y) origin 
       - oshape   - [in] ouitput image shape
    """

    sizex = arr.shape[0]
    sizey = arr.shape[1]

    #scale = float(L)

    xc, yc = center if center is not None else (sizex//2, sizey//2) 

    xarr = np.arange(math.floor(-xc), math.floor(sizex-xc))
    yarr = np.arange(math.floor(-yc), math.floor(sizey-yc))

    x,y = np.meshgrid(yarr, xarr) ### SWAPPED yarr, xarr to keep correct shape for grids
    r = np.sqrt(x*x+y*y)

    c = math.cos(math.radians(phi_deg))
    s = math.sin(math.radians(phi_deg))

    xrot = x*c - y*s
    yrot = x*s + y*c

    xrot = np.ceil(xrot)
    yrot = np.ceil(yrot)

    orows, orows1 = oshape[0], oshape[0] - 1
    ocols, ocols1 = oshape[1], oshape[1] - 1

    icols = np.array(xrot + math.ceil(ocols/2), dtype=np.int32)
    irows = np.array(yrot + math.ceil(orows/2), dtype=np.int32)

    irows = np.select([irows<0, irows>orows1], [0,orows1], default=irows)
    icols = np.select([icols<0, icols>ocols1], [0,ocols1], default=icols)

    sp.image = np.zeros(oshape, dtype=arr.dtype)
    sp.count = np.zeros(oshape, dtype=np.int32)

    unused_lst = list(map(_fillimg, irows, icols, arr))

    countpro = np.select([sp.count<1], [-1], default=sp.count)
    arr_rot = np.select([countpro>0], [sp.image/countpro], default=0)

    return arr_rot

#------------------------------

def plot_init() :
    import pyimgalgos.GlobalGraphics as gg
    global gg
    sp.fig, sp.axim, sp.axcb, sp.imsh = gg.fig_axim_axcb_imsh(figsize=(11,10))

#------------------------------

def plot_image(img) :

    ave, rms = img.mean(), img.std()
    amin, amax = ave-1*rms, ave+8*rms

    sp.axim.clear()
    if sp.imsh is not None : del sp.imsh
    sp.imsh = None

    gg.plot_imgcb(sp.fig, sp.axim, sp.axcb, sp.imsh, img, amin=amin, amax=amax, title=None)

    sp.fig.canvas.set_window_title('Test')
    sp.fig.canvas.draw() # re-draw figure content

    gg.show(mode='do not hold')

#------------------------------

def test_image_rotation() :
    from pyimgalgos.NDArrGenerators import random_standard

    #cent = (750, 750)
    shape = (1501,1501)
    img = random_standard(shape, mu=200, sigma=25, dtype=np.float32)

    plot_init()
    for phi in np.linspace(0,360,72, endpoint=False) :
        img_rot = image_rotation(img, phi, center=None, oshape=(2000,2000))
        print('test image rotation for phi(deg)=%.1f shape:%s' % (phi, str(img_rot.shape)))
        plot_image(img_rot)

    gg.show()

#------------------------------

if __name__ == "__main__" :
    test_image_rotation()
    sys.exit('The End')

#------------------------------
