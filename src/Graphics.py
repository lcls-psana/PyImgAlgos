#------------------------------
import numpy as np

import matplotlib
#if matplotlib.get_backend() != 'Qt4Agg' : matplotlib.use('Qt4Agg')

import matplotlib.pyplot  as plt
#import matplotlib.lines   as lines
#import matplotlib.patches as patches

from CalibManager.PlotImgSpeWidget import add_stat_text

#------------------------------
#class Storage :
#    def __init__(self) :
#        pass
#
#------------------------------
#store = Storage() # singleton
#------------------------------

#------------------------------

def figure(figsize=(13,12), title='Image', dpi=80, facecolor='w', edgecolor='w', frameon=True,
           move=None\
          ) :
    """ Creates and returns figure
    """
    fig = plt.figure(figsize=figsize,\
                     dpi=dpi,\
                     facecolor=facecolor,\
                     edgecolor=edgecolor,\
                     frameon=frameon)
    fig.canvas.set_window_title(title)
    if move is not None : move_fig(fig, x0=move[0], y0=move[1])
    return fig

#------------------------------

def move_fig(fig, x0=200, y0=100) :
    fig.canvas.manager.window.geometry('+%d+%d' % (x0, y0))

#------------------------------

def set_win_title(fig, titwin='Image') :
    fig.canvas.set_window_title(titwin)

#------------------------------

def add_axes(fig, axwin=(0.05, 0.03, 0.87, 0.93)) :
    """Add axes to figure from input list of windows.
    """
    axes = fig.add_axes(axwin)
    return axes

#------------------------------

def show(mode=None) :
    if mode is None : plt.ioff() # hold contraol at show() (connect to keyboard for controllable re-drawing)
    else            : plt.ion()  # do not hold control
    plt.show()

#------------------------------

def save_plt(fname='img.png', verb=True) :
    if verb : print 'Save plot in file: %s' % fname 
    plt.savefig(fname)

#------------------------------

def save_fig(fig, fname='img.png', verb=True) :
    if verb : print 'Save figure in file: %s' % fname 
    fig.savefig(fname)

#------------------------------

def hist(axhi, arr, bins=None, amp_range=None, weights=None, color=None, log=False) :
    """Makes historgam from input array of values (arr), which are sorted in number of bins (bins) in the range (amp_range=(amin,amax))
    """
    axhi.cla()
    hi = axhi.hist(arr.flatten(), bins=bins, range=amp_range, weights=weights, color=color, log=log) #, log=logYIsOn)
    if amp_range is not None : axhi.set_xlim(amp_range) # axhi.set_autoscale_on(False) # suppress autoscailing
    wei, bins, patches = hi
    add_stat_text(axhi, wei, bins)
    return hi

#------------------------------

def add_title_labels_to_axes(axes, title=None, xlabel=None, ylabel=None, fslab=14, fstit=20, color='k') :
    if title  is not None : axes.set_title(title, color=color, fontsize=fstit)
    if xlabel is not None : axes.set_xlabel(xlabel, fontsize=fslab)
    if ylabel is not None : axes.set_ylabel(ylabel, fontsize=fslab)

#------------------------------

def imshow_cbar(fig, axim, axcb, img, amin=None, amax=None, extent=None,\
                interpolation='nearest', aspect='auto', origin='upper',\
                orientation='horizontal') : # extent=img_range
    #if axim is None or axcb is None
    axim.cla()
    imsh = axim.imshow(img, interpolation=interpolation, aspect=aspect, origin=origin, extent=extent)
    cbar = fig.colorbar(imsh, cax=axcb, orientation=orientation)
    return imsh, cbar

#------------------------------

#------------------------------

#------------------------------

#------------------------------

#------------------------------

#------------------------------

#------------------------------

#------------------------------

