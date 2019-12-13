""" Loop over events in the dataset, get time in sec and nsec and feducial from EventId, make timestamp.
"""
from __future__ import print_function

import psana
from time import localtime, strftime

ds  = psana.DataSource('exp=cxif5315:run=169')
env = ds.env()

# loop over events in data set
for i, evt in enumerate(ds.events()) :

    if i>5 : break

    evtid = evt.get(psana.EventId)
    time_sec, time_nsec = evtid.time()
    tstamp = strftime('%Y-%m-%d %H:%M:%S', localtime(time_sec))
    print('ev#%2d   t=%d sec %d nsec    fid=%d    tstamp=%s' % (i, time_sec, time_nsec, evtid.fiducials(), tstamp))    
