#!/usr/bin/env python

##-----------------------------
from pyimgalgos.TDFileContainer import TDFileContainer
from pyimgalgos.TDIndexRecord   import TDIndexRecord
##-----------------------------

fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/lut-cxif5315-r0169-2015-10-28T15:32:20.txt'
fc = TDFileContainer(fname, indhdr='index', objtype=TDIndexRecord)#, pbits=1023)
fc.print_attrs()

for i, group in enumerate(fc.group_num_iterator()) :
    #if i>1000 : break
    group = fc.next()
    #group.print_attrs()    
    print '\n %s' % fc.hdr
    # Iterate over records in the group
    for rec in group() : # group() or group.get_objs()
        print rec.line.rstrip('\n')
        #rec.print_short()
        #print rec.index, rec.beta, rec.omega, rec.h, rec.k, rec.l, rec.dr, rec.R, rec.qv, rec.qh, rec.P

##-----------------------------
