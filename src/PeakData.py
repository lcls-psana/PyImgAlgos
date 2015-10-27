#------------------------------
"""Class helps to retreive and use peak data in processing

Usage::

    # Imports
    from pyimgalgos.PeakData import PeakData

    # Usage

    # make object
    pk = PeakData(line)

    # access peak attributes
    exp  = pk.exp  # (str)   experiment name
    run  = pk.run  # (int)   run number
    son  = pk.son  # (float) S/N for pixel with maximal intensity 
    sonc = pk.sonc # (float) S/N for all pixels included in the peak
    line = pk.line # (str)   entire record with peak data
    ...

    # get evaluated parameters
    signal = pk.peak_signal()
    noise  = pk.peak_noise()
    son    = pk.peak_son()

    # print attributes
    pk.print_peak_data()
    pk.print_peak_data_short()
    pk.print_attrs()

This software was developed for the SIT project.
If you use all or part of it, please give an appropriate acknowledgment.

@version $Id$

@author Mikhail S. Dubrovin
"""

#--------------------------------

import math
#import numpy as np
#from time import strftime, localtime #, gmtime

#------------------------------

class PeakData :

    def __init__(sp, line, pixel_size = 109.92) :  
        """Parse the string of parameters to values
        """
        ## Exp     Run  Date       Time      time(sec)   time(nsec) fiduc  Evnum  Reg  Seg  Row  Col  Npix      Amax      Atot   rcent   ccent rsigma  csigma rmin rmax cmin cmax    bkgd     rms     son  imrow   imcol     x[um]     y[um]     r[um]  phi[deg]
        #cxif5315  169  2015-02-22 02:20:47  1424600447  494719789  104424     1  EQU   17  170   51    38     168.6    2309.2   169.8    51.6   3.09    1.45  165  176   46   57   -2.90   27.99    6.12    586     499     -8027      -949      8082   -173.26

        sp.pixel_size = pixel_size
        
        sp.fields = line.rstrip('\n').split()

        s_exp, s_run, s_date, s_time, s_time_sec, s_time_nsec, \
        s_fid, s_evnum, s_reg, s_seg, s_row, s_col, s_npix, s_amax, s_atot, \
        s_rcent, s_ccent, s_rsigma, s_csigma, s_rmin, s_rmax, s_cmin, s_cmax, \
        s_bkgd, s_rms, s_son, s_imrow, s_imcol, s_x, s_y, s_r, s_phi = \
        sp.fields[0:32]
        
        sp.exp, sp.run, sp.evnum, sp.reg = s_exp, int(s_run), int(s_evnum), s_reg
        sp.date, sp.time, sp.tsec, sp.tnsec, sp.fid = s_date, s_time, int(s_time_sec), int(s_time_nsec), int(s_fid)
        sp.seg, sp.row, sp.col, sp.amax, sp.atot, sp.npix = int(s_seg), int(s_row), int(s_col), float(s_amax), float(s_atot), int(s_npix)
        sp.rcent, sp.ccent, sp.rsigma, sp.csigma = float(s_rcent), float(s_ccent), float(s_rsigma), float(s_csigma)
        sp.rmin, sp.rmax, sp.cmin, sp.cmax = int(s_rmin), int(s_rmax), int(s_cmin), int(s_cmax)
        sp.bkgd, sp.rms, sp.son = float(s_bkgd), float(s_rms), float(s_son)
        sp.imrow, sp.imcol = int(s_imrow), int(s_imcol)
        sp.x, sp.y, sp.r, sp.phi = float(s_x), float(s_y), float(s_r)/sp.pixel_size, float(s_phi)
        sp.sonc = sp.peak_son()
        sp.dphi000 = sp.phi
        sp.dphi180 = sp.phi - 180 if sp.phi>-90 else sp.phi + 180 # +360-180

        sp.line = line
        
#------------------------------

    def print_peak_data_short(sp) :
        """Prints short subset of data
        """    
        print '%5d %s %3d %3d %3d %7.1f %7.1f %3d %6d %6d %7.1f %7.1f' % \
              (sp.evnum, sp.reg, sp.seg, sp.row, sp.col, sp.amax, sp.atot, sp.npix, sp.x, sp.y, sp.r, sp.phi)   

#------------------------------

    def print_peak_data(sp) :
        """Prints input data string(line)
        """    
        for field in sp.fields : print field,
        print ''

#------------------------------

    def peak_signal(sp) :
        """Evaluates corrected signal subtracting the background
        """
        return sp.atot-sp.bkgd*sp.npix

#------------------------------

    def peak_noise(sp) :
        """Evaluates corrected rms noise for all pixels in the peak
        """
        return sp.rms*math.sqrt(sp.npix)

#------------------------------

    def peak_son(sp) :
        """Evaluates corrected value of the S/N ratio based on entire peak intensity
        """
        N = sp.peak_noise()
        return sp.peak_signal()/N if N>0 else 0

#------------------------------

    def print_attrs(sp) :
        msg = 'Attributes of %s, pixel size[um] =%8.2f' % (sp.__class__.__name__, sp.pixel_size)
        #msg += ', line:  \n%s' % (sp.line)
        print msg

#------------------------------
#--------------------------------
#-----------  TEST  -------------
#--------------------------------

if __name__ == "__main__" :
    pass

#--------------------------------
