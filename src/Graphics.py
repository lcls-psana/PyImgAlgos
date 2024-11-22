#!@PYTHON@
####!/usr/bin/env python

""":py:class:`Graphics` wrapping methods for matplotlib
=======================================================

Usage::

    import pyimgalgos.Graphics as gr

    # Methods

    fig = gr.figure(figsize=(13,12), title='Image', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=None)
    gr.move_fig(fig, x0=200, y0=100)
    gr.move(x0=200, y0=100)
    gr.add_axes(fig, axwin=(0.05, 0.03, 0.87, 0.93))
    gr.fig_img_cbar_axes(fig=None, win_axim=(0.05,  0.03, 0.87, 0.93), win_axcb=(0.923, 0.03, 0.02, 0.93))
    gr.set_win_title(fig, titwin='Image')
    gr.add_title_labels_to_axes(axes, title=None, xlabel=None, ylabel=None, fslab=14, fstit=20, color='k')
    gr.show(mode=None)
    gr.draw()
    gr.draw_fig(fig)
    gr.save_plt(fname='img.png', verb=True)
    gr.save_fig(fig, fname='img.png', verb=True)
    hi = gr.hist(axhi, arr, bins=None, amp_range=None, weights=None, color=None, log=False)
    imsh = gr.imshow(axim, img, amp_range=None, extent=None, interpolation='nearest', aspect='auto', origin='upper', orientation='horizontal', cmap='inferno')
    cbar = gr.colorbar(fig, imsh, axcb, orientation='vertical', amp_range=None)
    imsh, cbar = gr.imshow_cbar(fig, axim, axcb, img, amin=None, amax=None, extent=None, interpolation='nearest', aspect='auto', origin='upper', orientation='vertical', cmap='inferno')

See:
  - :py:class:`Graphics`
  - :py:class:`GlobalGraphics`
  - :py:class:`NDArrGenerators`
  - :py:class:`HBins`
  - :py:class:`HPolar`
  - :py:class:`HSpectrum`
  - :py:class:`NDArrSpectrum`
  - :py:class:`RadialBkgd`
  - `Radial background <https://confluence.slac.stanford.edu/display/PSDMInternal/Radial+background+subtraction+algorithm>`_.
  - `matplotlib <https://matplotlib.org/contents.html>`_.

This software was developed for the SIT project.
If you use all or part of it, please give an appropriate acknowledgment.

Created in 2015 by Mikhail Dubrovin
"""
from __future__ import print_function

import numpy as np

import matplotlib
#if matplotlib.get_backend() != 'Qt4Agg': matplotlib.use('Qt4Agg')

import matplotlib.pyplot  as plt
#import matplotlib.lines   as lines
#import matplotlib.patches as patches

from CalibManager.PlotImgSpeWidget import add_stat_text


def dict_subset(d, keys):
    return {k:v for k,v in d.items() if k in keys}


def set_win_title(fig, titwin='Image'):
    fig.canvas.manager.set_window_title(titwin)


def figure(**kwa):
    """ Creates and returns figure.
        figsize=(13,12), title='Image', dpi=80, facecolor='w', edgecolor='w', frameon=True
    """
    fig = plt.figure(**dict_subset(kwa, ('num', 'figsize', 'dpi', 'facecolor', 'edgecolor', 'frameon', 'FigureClass', 'clear',\
                                         'linewidth', 'subplotpars', 'tight_layout', 'constrained_layout')))
    move = kwa.get('move', None)
    title = kwa.get('title', '')
    if title: set_win_title(fig, title)
    if move: move_fig(fig, x0=move[0], y0=move[1])
    return fig


def pp_hist(axis, x, **kwa):
    """ matplotlib.pyplot.hist(x,
                       bins=10,
                       range=None,
                       normed=False,
                       weights=None,
                       cumulative=False,
                       bottom=None,
                       histtype=u'bar',
                       align=u'mid',
                       orientation=u'vertical',
                       rwidth=None,
                       log=False,
                       color=None,
                       label=None,
                       stacked=False,
                       hold=None,
                       **kwargs)
    """
    return axis.hist(x, **dict_subset(kwa,\
           ('bins', 'range', 'normed', 'weights', 'cumulative', 'bottom', 'histtype', 'align',\
            'orientation', 'rwidth', 'log', 'color', 'label', 'stacked', 'hold')))


def move_fig(fig, x0=200, y0=100):
    #fig.canvas.manager.window.geometry('+%d+%d' % (x0, y0)) # in previous version of matplotlib
    backend = matplotlib.get_backend()
    #logger.debug('matplotlib.get_backend(): %s' % backend)
    if backend == 'TkAgg': # this is our case
        fig.canvas.manager.window.wm_geometry("+%d+%d" % (x0, y0))
    elif backend == 'WXAgg':
        fig.canvas.manager.window.SetPosition((x0, y0))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        fig.canvas.manager.window.move(x0, y0)


def move(x0=200,y0=100):
    move_fig(plt.gcf(), x0, y0)
    #plt.get_current_fig_manager().window.move(x0, y0)
    #plt.get_current_fig_manager().window.geometry('+%d+%d' % (x0, y0))


def add_axes(fig, axwin=(0.05, 0.03, 0.87, 0.93)):
    """Add axes to figure from input list of windows.
    """
    return fig.add_axes(axwin)


def fig_img_axes(fig=None, win_axim=(0.08,  0.05, 0.89, 0.93)):
    """ Returns figure and image axes
    """
    _fig = figure(figsize=(6,5)) if fig is None else fig
    axim = _fig.add_axes(win_axim)
    return _fig, axim


def fig_axes(fig, windows=((0.05,  0.03, 0.87, 0.93), (0.923, 0.03, 0.02, 0.93))):
    """ Returns list of figure axes for input list of windows
    """
    return [fig.add_axes(w) for w in windows]


def fig_img_cbar_axes(fig=None,\
             win_axim=(0.05,  0.03, 0.87, 0.93),\
             win_axcb=(0.923, 0.03, 0.02, 0.93)):
    """ Returns figure and axes for image and color bar
    """
    _fig = figure() if fig is None else fig
    axim = _fig.add_axes(win_axim)
    axcb = _fig.add_axes(win_axcb)
    return _fig, axim, axcb


def fig_axim_axcb_imsh(figsize=(13,12), title='Image', dpi=80,\
                       win_axim=(0.05,  0.03, 0.87, 0.93),\
                       win_axcb=(0.923, 0.03, 0.02, 0.93),\
                       arr2d=np.zeros((10,10)), origin='upper'):
    """ Creates and returns figure, axes for image and color bar, imshow object
    """
    fig  = plt.figure(figsize=figsize, dpi=dpi, facecolor='w', edgecolor='w', frameon=True)
    axim = fig.add_axes(win_axim)
    axcb = fig.add_axes(win_axcb)
    set_win_title(fig, title)
    imsh = axim.imshow(arr2d, interpolation='nearest', aspect='auto', origin=origin)
    return fig, axim, axcb, imsh


FYMIN, FYMAX = 0.050, 0.90
def fig_img_cbar_hist_axes(fig=None,\
                      win_axim=(0.02,  FYMIN, 0.8,  FYMAX),\
                      win_axcb=(0.915, FYMIN, 0.01, FYMAX),\
                      win_axhi=(0.76,  FYMIN, 0.15, FYMAX),\
                      **kwa):
    """ Returns figure and axes for image, color bar, and spectral histogram
    """
    _fig = figure() if fig is None else fig
    return _fig,\
           _fig.add_axes(win_axim, **kwa),\
           _fig.add_axes(win_axcb, **kwa),\
           _fig.add_axes(win_axhi, **kwa)


def add_title_labels_to_axes(axes, title=None, xlabel=None, ylabel=None, fslab=14, fstit=20, color='k'):
    if title  is not None: axes.set_title(title, color=color, fontsize=fstit)
    if xlabel is not None: axes.set_xlabel(xlabel, fontsize=fslab)
    if ylabel is not None: axes.set_ylabel(ylabel, fontsize=fslab)


def show(mode=None):
    #plt.hold(True)
    if mode is None: plt.ioff() # hold contraol at show() (connect to keyboard for controllable re-drawing)
    else           : plt.ion()  # do not hold control
    plt.pause(0.001) # hack to make it work... othervise show() does not work...
    plt.show()


def draw():
    plt.draw()


def draw_fig(fig):
    fig.canvas.draw()


def save_plt(fname='img.png', verb=True):
    if verb: print('Save plot in file: %s' % fname)
    plt.savefig(fname)


def save_fig(fig, fname='img.png', verb=True):
    if verb: print('Save figure in file: %s' % fname)
    fig.savefig(fname)


def hist(axhi, arr, bins=None, amp_range=None, weights=None, color=None, log=False):
    """Makes historgam from input array of values (arr), which are sorted in number of bins (bins) in the range (amp_range=(amin,amax))
    """
    #axhi.cla()
    hi = axhi.hist(arr.flatten(), bins=bins, range=amp_range, weights=weights, color=color, log=log) #, log=logYIsOn)
    if amp_range is not None: axhi.set_xlim(amp_range) # axhi.set_autoscale_on(False) # suppress autoscailing
    wei, bins, patches = hi
    add_stat_text(axhi, wei, bins)
    return hi


def imshow(axim, img, amp_range=None, extent=None,\
           interpolation='nearest', aspect='auto', origin='upper',\
           orientation='horizontal', cmap='inferno'):
    """
    extent - list of four image physical limits for labeling,
    cmap: 'jet', 'gray_r', 'inferno'
    #axim.cla()
    """
    imsh = axim.imshow(img, interpolation=interpolation, aspect=aspect, origin=origin, extent=extent, cmap=cmap)
    if amp_range is not None: imsh.set_clim(amp_range[0],amp_range[1])
    return imsh


def colorbar(fig, imsh, axcb, orientation='vertical', amp_range=None):
    """
    orientation = 'horizontal'
    amp_range = (-10,50)
    """
    if amp_range is not None: imsh.set_clim(amp_range[0],amp_range[1])
    cbar = fig.colorbar(imsh, cax=axcb, orientation=orientation)
    return cbar


def imshow_cbar(fig, axim, axcb, img, amin=None, amax=None, extent=None,\
                interpolation='nearest', aspect='auto', origin='upper',\
                orientation='vertical', cmap='inferno'):
    """
    extent - list of four image physical limits for labeling,
    cmap: 'gray_r'
    #axim.cla()
    """
    axim.cla()
    if img is None: return
    ave = np.mean(img) if amin is None and amax is None else None
    rms = np.std(img)  if amin is None and amax is None else None
    cmin = amin if amin is not None else ave-1*rms if ave is not None else None
    cmax = amax if amax is not None else ave+3*rms if ave is not None else None

    imsh = axim.imshow(img, interpolation=interpolation, aspect=aspect, origin=origin, extent=extent, cmap=cmap, vmin=cmin, vmax=cmax)
    cbar = fig.colorbar(imsh, cax=axcb, orientation=orientation)
    if cmin is not None:
      imsh.set_clim(cmin, cmax)
      #cbar.set_clim(cmin, cmax)
    return imsh, cbar


def test01():
    """ imshow
    """
    img = random_standard(shape=(40,60), mu=200, sigma=25)
    #fig = figure(figsize=(6,5), title='Test imshow', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=(100,10))
    #axim = add_axes(fig, axwin=(0.10, 0.08, 0.85, 0.88))
    fig, axim = fig_img_axes()
    move_fig(fig, x0=200, y0=100)
    imsh = imshow(axim, img, amp_range=None, extent=None,\
           interpolation='nearest', aspect='auto', origin='upper',\
           orientation='horizontal', cmap='jet')


def test02():
    """ hist
    """
    mu, sigma = 200, 25
    arr = random_standard((500,), mu, sigma)
    #fig = figure(figsize=(6,5), title='Test hist', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=(100,10))
    #axhi = add_axes(fig, axwin=(0.10, 0.08, 0.85, 0.88))
    fig, axhi = fig_img_axes()
    his = hist(axhi, arr, bins=100, amp_range=(mu-6*sigma,mu+6*sigma), weights=None, color=None, log=False)


def test03():
    """ Update image in the event loop
    """
    #fig = figure(figsize=(6,5), title='Test hist', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=(100,10))
    #axim = add_axes(fig, axwin=(0.10, 0.08, 0.85, 0.88))
    fig, axim = fig_img_axes()
    imsh = None
    for i in range(10):
       print('Event %3d' % i)
       img = random_standard((1000,1000), mu=200, sigma=25)
       #axim.cla()
       set_win_title(fig, 'Event %d' % i)

       if imsh is None:
           imsh = imshow(axim, img, amp_range=None, extent=None,\
                  interpolation='nearest', aspect='auto', origin='upper',\
                  orientation='horizontal', cmap='jet')
       else:
           imsh.set_data(img)
       show(mode=1)
       #draw_fig(fig)


def test04():
    """ Update histogram in the event loop
    """
    mu, sigma = 200, 25
    #fig = figure(figsize=(6,5), title='Test hist', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=(100,10))
    #axhi = add_axes(fig, axwin=(0.10, 0.08, 0.85, 0.88))
    fig, axhi = fig_img_axes()

    for i in range(10):
       print('Event %3d' % i)
       arr = random_standard((500,), mu, sigma, dtype=np.float32)
       axhi.cla()
       set_win_title(fig, 'Event %d' % i)
       his = hist(axhi, arr, bins=100, amp_range=(mu-6*sigma,mu+6*sigma), weights=None, color=None, log=False)

       show(mode=1)
       #draw(fig)


def test05():
    """ Update image with color bar in the event loop
    """
    fig, axim, axcb = fig_img_cbar_axes()
    move_fig(fig, x0=200, y0=0)
    imsh = None
    for i in range(20):
       print('Event %3d' % i)
       img = random_standard((1000,1000), mu=i, sigma=10)
       #axim.cla()
       set_win_title(fig, 'Event %d' % i)
       if imsh is None:
           imsh, cbar = imshow_cbar(fig, axim, axcb, img, amin=None, amax=None, extent=None,\
                                    interpolation='nearest', aspect='auto', origin='upper',\
                                    orientation='vertical', cmap='inferno')
       else:
           imsh.set_data(img)
           ave, rms = img.mean(), img.std()
           imsh.set_clim(ave-1*rms, ave+3*rms)
       show(mode=1)
       #draw_fig(fig)


def test_selected():

    from time import time
    import sys; global sys
    from pyimgalgos.NDArrGenerators import random_standard; global random_standard

    if len(sys.argv)==1:
        print('Use command > python %s <test-number [1-5]>' % sys.argv[0])
        sys.exit ('Add <test-number> in command line...')

    tname = sys.argv[1] if len(sys.argv) > 1 else '1'
    print(50*'_', '\nTest %s' % tname)

    t0_sec=time()
    if   tname == '1': test01()
    elif tname == '2': test02()
    elif tname == '3': test03()
    elif tname == '4': test04()
    elif tname == '5': test05()
    else: sys.exit('Test %s is not implemented' % tname)
    msg = 'Test %s consumed time %.3f' % (tname, time()-t0_sec)
    show()
    sys.exit(msg)


def test_all():
    test01()
    test02()
    show()


if __name__ == "__main__":
    test_selected()
    #test_all()
    print('End of test')

# EOF

