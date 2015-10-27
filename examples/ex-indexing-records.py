
# Import
from pyimgalgos.TDFileContainer import TDFileContainer
from pyimgalgos.TDIndexRecord  import TDIndexRecord

# Initialization
# fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/lut-cxif5315-r0169-2015-10-23T16:23:44.txt'
fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/lut-cxif5315-r0169-2015-10-26T16:56:41.txt'
fc = TDFileContainer(fname, indhdr='index', objtype=TDIndexRecord)#, pbits=1023)
fc.print_attrs()

#print 'lst_grnum   : ', fc.lst_grnum
#print 'lst_begin   : ', fc.lst_begin
#print 'lst_nrecords: ', fc.lst_nrecords

# Iterate over groups
for i, group in enumerate(fc.group_num_iterator()) :
    group = fc.next()
    #group.print_attrs()
    lst_recs = group() # group() or group.get_objs()
    
    #if i>1000 : break

    #if len(lst_recs) :
    print '\n %s' % fc.hdr
    # Iterate over records in the group
    for rec in lst_recs :
        rec.print_short()
