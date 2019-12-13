#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import pyimgalgos.GlobalGraphics as gg

#------------------------------

def do_work(path='img-averaged.npy') :
    print('Load image from file %s' % path)
    img = np.load(path)
    #ave, rms = img.mean(), img.std()
    #amin, amax = ave-1*rms, ave+10*rms
    amin, amax = -10, 80
    fig, axim, axcb = gg.fig_axes(figsize=(13,12), title='Image', dpi=80, win_axim=(0.05,  0.03, 0.87, 0.93), win_axcb=(0.923, 0.03, 0.02, 0.93)) 
    gg.plot_img(img, amin=amin, amax=amax)
    #plotImageLarge(arr, img_range=None, amp_range=None, figsize=(12,10), title='Image', origin='upper', window=(0.05,  0.03, 0.94, 0.94))
    gg.show()

    ofname = '%s.png' % os.path.splitext(path)[0]
    gg.save_fig(fig, ofname, pbits=1)

#------------------------------

def usage() :
    print('Usage: %s <image-file-name>.npy' % str(sys.argv[0]))
    
#------------------------------

if __name__ == '__main__' :
    usage()
    if len(sys.argv)>1 : do_work(sys.argv[1])
    else               : do_work()
    sys.exit ('End of script')

#------------------------------
