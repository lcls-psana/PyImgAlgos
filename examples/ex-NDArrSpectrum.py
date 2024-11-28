#!/usr/bin/env python

from time import time, sleep
import psana
import numpy as np
from Detector.AreaDetector import AreaDetector
from Detector.GlobalUtils import print_ndarr
import pyimgalgos.Graphics as gr
from pyimgalgos.NDArrSpectrum import NDArrSpectrum

DO_PLOT = False
#DO_PLOT = True

runnum = 54
ds  = psana.DataSource('exp=xpptut15:run=%d' % runnum)
env = ds.env()
src = psana.Source('XppGon.0:Cspad.0')

det = AreaDetector(src, env, pbits=0)
#evt = ds.events().next()

peds  = det.pedestals(runnum)
smask = det.status_as_mask(runnum)

range = (-20, 80)
nbins = 100
spec = NDArrSpectrum(range, nbins)

if DO_PLOT :
    fighi = gr.figure(figsize=(6,5), title='Spectrum')
    axhi  = gr.add_axes(fighi, axwin=(0.15, 0.10, 0.82, 0.85))
    gr.move_fig(fighi, x0=800, y0=0)

    figim = gr.figure(figsize=(13,12), title='Image')
    axim  = gr.add_axes(figim, axwin=(0.05, 0.03, 0.87, 0.93))
    gr.move_fig(figim, x0=0, y0=0)

arr_dt = []

for i,evt in enumerate(ds.events()) :
    if i>10 : break

    t0_sec = time()
    raw = det.raw(evt) # dtype:int16
    # print_ndarr(raw, 'raw')
    nda = np.array(raw, dtype=np.float32)
    if nda is None : continue
    nda -= peds
    det.common_mode_apply(runnum, nda)
    nda *= smask

    spec.fill(nda)

    dt = time()-t0_sec
    print('Event %3d processing time %.3f sec' % (i, dt))
    arr_dt.append(dt)

    if DO_PLOT :

        img = det.image(evt,nda)

        hi = gr.hist(axhi, nda, bins=100, amp_range=(-50, 50))
        gr.set_win_title(fighi, titwin='Spectrum for event:%3d'%i)

        ave, rms = nda.mean(), nda.std()
        imsh = axim.imshow(img) #, interpolation=interpolation, aspect=aspect, origin=origin, extent=extent)
        imsh.set_clim(ave-1*rms, ave+2*rms)
        gr.set_win_title(figim, titwin='Image for event:%3d'%i)

        fighi.canvas.draw()
        figim.canvas.draw()
        gr.show(mode='go')

if DO_PLOT : gr.show(mode=None)

histarr, edges, nbins = spec.spectrum()
print_ndarr(histarr, 'histarr')

nparr_dt = np.array(arr_dt)

print('Number of events with data processed = %d' % len(arr_dt))
print('Average time per event dt/event=%.3f sec, RMS(dt)=%.3f sec' % (nparr_dt.mean(), nparr_dt.std()))

#------------------------------
