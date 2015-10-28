#------------------------------
"""Class helps to retreive and use peak data in processing

Usage::

    # Imports
    from pyimgalgos.TDIndexRecord import TDIndexRecord

    # Usage

    # make object
    rec = TDIndexRecord(line)

    # access record attributes
    index, beta, omega h, k, l, dr, R, qv, qh, P =\
    rec.index, rec.beta, rec.omega rec.h, rec.k, rec.l, rec.dr, rec.R, rec.qv, rec.qh, rec.P
    line = rec.line

    # print attributes
    rec.print_short()

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

class TDIndexRecord :

    def __init__(sp, line) : # , sigma_q = 0.000484) :  
        """Parse the string of parameters to values

        # sigma_q   = 0.000484 1/A (approximately pixel size/sample-to-detector distance = 100um/100mm)

        # beta 10.00  omega 0.50 degree
        # index   beta     omega   h  k  l     dr [1/A]   R(h,k,l)   qv [1/A]   qh [1/A]   P(omega)
             2    10.00     0.50   0  1  0     0.001197   0.038484  -0.000058   0.038418   0.047045
        """
        sp.fields = line.rstrip('\n').split()

        s_index, s_beta, s_omega, s_h, s_k, s_l, s_dr, s_R, s_qv, s_qh, s_P = sp.fields[0:11]
        sp.index, sp.beta, sp.omega = int(s_index), float(s_beta), float(s_omega)
        sp.h, sp.k, sp.l = int(s_h), int(s_k), int(s_l)
        sp.dr, sp.R, sp.qv, sp.qh, sp.P = float(s_dr), float(s_R), float(s_qv), float(s_qh), float(s_P) 

        sp.line = line
        
#------------------------------

    def print_short(sp) :
        """Prints short subset of data
        """    
        print '%6d  %7.2f  %7.2f  %2d %2d %2d    %9.6f  %9.6f  %9.6f  %9.6f  %9.6f' % \
              (sp.index, sp.beta, sp.omega, sp.h, sp.k, sp.l, sp.dr, sp.R, sp.qv, sp.qh, sp.P)   

#------------------------------
#--------------------------------
#-----------  TEST  -------------
#--------------------------------

if __name__ == "__main__" :
    pass

#--------------------------------
