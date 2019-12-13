#!/usr/bin/env python

##-----------------------------
from __future__ import print_function
from pyimgalgos.TDFileContainer import TDFileContainer
from time import time
#from pyimgalgos.TDPeak          import TDPeak - is used by default
##-----------------------------

def do_work() :

    fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/pfv2-cxif5315-r0169-2015-09-14T14:28:04.txt'
    fc = TDFileContainer(fname, indhdr='Evnum') #, objtype=TDPeak, pbits=0)
    #fc.print_attrs()
    fc.print_content(nlines=20) # print part of content 

    # Direct access to TDGroup object
    group = fc.group(8)
    group.print_attrs()

    t0_sec = time()

    for grpnum in fc.group_num_iterator() :
        group = next(fc)
        print('%s\nEvent %d  ' % (80*'_', grpnum))
        for i, pk in enumerate(group()) :
            print('  peak#%2d  %s  %s  t[nsec]=%8d  seg=%2d  row=%3d  col=%3d  %s  bkgd=%5.1f  rms=%5.1f  S/N=%5.1f  r=%6.1f  phi=%6.1f' % \
                  (i, pk.exp, pk.time, pk.tnsec, pk.seg, pk.row, pk.col, pk.reg, pk.bkgd, pk.rms, pk.son, pk.r, pk.phi))
            #rec.print_short()

            # Information available through the PeakData object pk
            # ____________________________________________________
            # pk.exp, pk.run, pk.evnum, pk.reg
            # pk.date, pk.time, pk.tsec, pk.tnsec, pk.fid
            # pk.seg, pk.row, pk.col, pk.amax, pk.atot, pk.npix
            # pk.rcent, pk.ccent, pk.rsigma, pk.csigma
            # pk.rmin, pk.rmax, pk.cmin, pk.cmax
            # pk.bkgd, pk.rms, pk.son
            # pk.imrow, pk.imcol
            # pk.x, pk.y, pk.r, pk.phi
            # pk.sonc
            # pk.dphi000
            # pk.dphi180
            # pk.line

            # get evaluated parameters
            # pk.peak_signal()
            # pk.peak_noise()
            # pk.peak_son()

            # print attributes
            # pk.print_peak_data()
            # pk.print_peak_data_short()
            # pk.print_attrs()
        
    print('\nTime to iterate using next() %.3f sec' % (time()-t0_sec))

    #t0_sec = time()
    #groups = fc.list_of_groups()  
    #print 'Time to generate list of group objects %.3f sec' % (time()-t0_sec)

##-----------------------------
if __name__ == "__main__" :
    do_work()
    print('Test is completed')
    #sys.exit('Processing is completed')
##-----------------------------
