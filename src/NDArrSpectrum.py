#!/usr/bin/env python

#--------------------------------
""":py:class:`NDArrSpectrum` - ALIAS FOR HSpectrum - support creation of spectral histogram for arbitrary shaped numpy array.

Usage::

    # Import
    # ==============
    from pyimgalgos.NDArrSpectrum import NDArrSpectrum


    # Initialization
    # ==============
    # 1) for bins of equal size:
    range = (vmin, vmax)
    nbins = 100
    spec = NDArrSpectrum(range, nbins)

    # 2) for variable size bins:
    bins = (v0, v1, v2, v4, v5, vN) # any number of bin edges
    spec = NDArrSpectrum(bins)


    # Fill spectrum
    # ==============
    # nda = ... (get it for each event somehow)
    spec.fill(nda)


    # Get spectrum
    # ==============
    histarr, edges, nbins = spec.spectrum()


    # Optional
    # ==============
    spec.print_attrs()

@see :py:class:`pyimgalgos.HSpectrum`,
:py:class:`pyimgalgos.HBins`
:py:class:`pyimgalgos.HPolar`

This software was developed for the SIT project.
If you use all or part of it, please give an appropriate acknowledgment.

Revision: $Revision$

@version $Id$

@author Mikhail S. Dubrovin

"""
#--------------------------------
__version__ = "$Revision$"
#--------------------------------

from pyimgalgos.HSpectrum import HSpectrum

NDArrSpectrum = HSpectrum

#------------------------------

if __name__ == "__main__" :
    from pyimgalgos.HSpectrum import main
    main()

#------------------------------
