#!/usr/bin/env python

##-----------------------------
from __future__ import print_function
from pyimgalgos.TDFileContainer import TDFileContainer
from pyimgalgos.TDIndexRecord   import TDIndexRecord
##-----------------------------

fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/lut-cxif5315-r0169-2015-10-28T15:32:20.txt'
fc = TDFileContainer(fname, indhdr='index', objtype=TDIndexRecord)#, pbits=1023)
fc.print_attrs()
fc.print_content(nlines=50)

for i, grpnum in enumerate(fc.group_num_iterator()) :
    #if i>1000 : break
    group = next(fc)
    #group.print_attrs()    
    print('%s\nOrientation group %d\n %s' % (80*'_', grpnum, fc.hdr))
    # Iterate over records in the group
    for rec in group() : # group() or group.get_objs()
        print(rec.line.rstrip('\n'))
        #rec.print_short()
        #print rec.index, rec.beta, rec.omega, rec.h, rec.k, rec.l, rec.dr, rec.R, rec.qv, rec.qh, rec.P

##-----------------------------
