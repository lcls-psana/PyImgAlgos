
"""
Author: Chris O'Grady
"""
from __future__ import print_function

from psana import *
import numpy
import time

class SmallData(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def getPsanaEvent(cheetahSmallData):
    hrsMinSec = cheetahSmallData.Filename.split('_')[4]
    fid = int(cheetahSmallData.Filename.split('_')[5].split('.')[0],16)
    for t in times:
        if t.fiducial() == fid:
            localtime = time.strftime('%H:%M:%S', time.localtime(t.seconds()))
            localtime = localtime.replace(':','')
            if localtime[0:3] == hrsMinSec[0:3]: return run.event(t)
    return None

def getSmallData(file):
    f = open(file)
    smldata=[]
    cheetahFile=True
    for line in f:
        if line.startswith('#'):
            if line.count(',') == 0: cheetahFile=False
            if cheetahFile:
                print('*** Found cheetah small data file')
            else:
                print('*** Found psana small data file')
            line = line[2:-1]
            line = line.replace('eventData->','')
            fieldNames = line.split()
            # hack some of the fieldNames
            removeChars = [',','[',']','(',')']
            for rc in removeChars:
                fieldNames = [fname.replace(rc,'') for fname in fieldNames]
            print('Small data fieldnames:',fieldNames)
        else:
            fields = line.split()
            fields = [field.replace(',','') for field in fields]
            event = SmallData()
            for i,field in enumerate(fields):
                try:
                    field = int(field)
                except ValueError:
                    try:
                        field = float(field)
                    except ValueError:
                        pass
                setattr(event,fieldNames[i],field)
            smldata.append(event)
    f.close()
    return cheetahFile,smldata

######################################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("exprun", help="psana experiment/run string (e.g. exp=xppd7114:run=43)")
parser.add_argument("smldatafile", help="cheetah or psana text file with hit information")
parser.add_argument("-c","--cfg",help="psana configuration file")
args = parser.parse_args()

cheetahFile,smallData = getSmallData(args.smldatafile)

if args.cfg:
    setConfigFile(args.cfg)

ds = DataSource('%s:idx'%args.exprun)
run = next(ds.runs())
times = run.times()

for sd in smallData:
    if cheetahFile:
        if sd.frameNumber==5:
            evt = getPsanaEvent(sd)
            if evt is None:
                print('*** Failed to find event')
            else:
                print('*** Jumped to analyze event with ID:',evt.get(EventId))
            break
    else:
        if sd.Npix==56:
            et = EventTime((sd.timesec<<32)|sd.timensec,sd.fiduc)
            evt = run.event(et)
            if evt is None:
                print('*** Failed to find event')
            else:
                print('*** Jumped to analyze event with ID:',evt.get(EventId))
            break
