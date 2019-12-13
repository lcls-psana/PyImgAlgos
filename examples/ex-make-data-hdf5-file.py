from __future__ import print_function

#------------------------------
# https://confluence.slac.stanford.edu/display/PSDM/MPI+Parallelization
#------------------------------
import psana
#------------------------------

def make_data_hdf5_t(tname='0', dsname='exp=xpptut15:run=54:smd', ofname='run54.h5') :

    dsource  = psana.MPIDataSource(dsname)
    cspaddet = psana.Detector('cspad')
    smldata  = dsource.small_data(ofname, gather_interval=100)
    run = next(dsource.runs())

    data_sum = None
    for nevt,t in enumerate(run.times()):
       evt = run.event(t)
       data = cspaddet.raw(evt) if tname in ['0','1']\
       else   cspaddet.calib(evt)

       if data is None: continue
 
       # save per-event data
       smldata.event(cspad_data=data)

       if data_sum is None: data_sum = data
       else:                data_sum += data
 
       if nevt>9: break
 
    # get "summary" data
    run_sum = smldata.sum(data_sum)
    # save HDF5 file, including summary data
    smldata.save(run_sum=run_sum)

#------------------------------

def make_data_hdf5(tname='0', dsname='exp=xpptut15:run=54:smd', ofname='run54.h5', det='CxiDs1.0:Cspad.0') :

    dsource  = psana.MPIDataSource(dsname)
    cspaddet = psana.Detector(det)
    smldata  = dsource.small_data(ofname, gather_interval=100)

    data_sum = None
    for nevt,evt in enumerate(dsource.events()):

        if nevt%1000==0 : print('Event %d'%nevt)

        if tname == '3' : # for amo86615:run=197 (Chuck) 
          if nevt>78141: break
          if not(nevt in (264, 58452, 78140, 76760)) : continue

        elif tname == '4' :  # for cxitut13:run=10 (Chuck)
          if nevt>212: break
          if not(nevt in (18, 50, 78, 211)) : continue

        else :
          if nevt>9: break


        data = cspaddet.raw(evt) if tname in ['0','1']\
        else   cspaddet.calib(evt).astype(np.int16)

        if data is None: continue

        # save per-event data
        smldata.event(cspad_data=data)

        if data_sum is None: data_sum  = data
        else:                data_sum += data
 
 
    # get "summary" data
    run_sum = smldata.sum(data_sum)
    # save HDF5 file, including summary data
    smldata.save(run_sum=run_sum)

#------------------------------

if __name__ == "__main__" :
    import sys; global sys
    import numpy as np; global np
    tname = sys.argv[1] if len(sys.argv) > 1 else '0'
    print(50*'_', '\nTest %s' % tname)
    if   tname == '0': make_data_hdf5(tname, 'exp=cxif5315:run=129:smd', 'smd-cxif5315-r129-dark.h5',        det='CxiDs2.0:Cspad.0')
    elif tname == '1': make_data_hdf5(tname, 'exp=cxif5315:run=169:smd', 'smd-cxif5315-r169-raw.h5',         det='CxiDs2.0:Cspad.0')
    elif tname == '2': make_data_hdf5(tname, 'exp=cxif5315:run=169:smd', 'smd-cxif5315-r169-calib-fde.h5',   det='CxiDs2.0:Cspad.0')
    elif tname == '3': make_data_hdf5(tname, 'exp=amo86615:run=197:smd', 'smd-amo86615-r197-calib-spi.h5',   det='Camp.0:pnCCD.0')
    elif tname == '4': make_data_hdf5(tname, 'exp=cxitut13:run=10:smd',  'smd-cxitut13-r010-calib-cryst.h5', det='CxiDs1.0:Cspad.0')
    else : sys.exit('Test %s is not implemented' % tname)
    sys.exit('End of Test %s' % tname)

#------------------------------
