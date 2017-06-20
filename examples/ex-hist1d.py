#!/usr/bin/env python

from pyimgalgos.GlobalGraphics import hist1d, show
import numpy as np

arr = 10 + 20*np.random.standard_normal((1000,))

# Returns plt.figure, fig.add_axes(axwin), axhi.hist
fig, axhi, hi =\
hist1d(arr, bins=None, amp_range=None, weights=None, color='y', show_stat=True, log=False,\
       figsize=(6,5), axwin=(0.15, 0.12, 0.78, 0.80),\
       title='My Title', xlabel='my x-label', ylabel='my x-label', titwin='ex-hist1d')

show()

