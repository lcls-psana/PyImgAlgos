#!/usr/bin/env python
"""`HBins.py` histogram-style bin parameters holder

Usage::

    from pyimgalgos.HBins import HBins

    # Equal bins
    hb = HBins((1,6), nbins=5)

    # Variable bins
    hb = HBins((1,2,4,6,10))

    # Access methods
    nbins         = hb.nbins()         # int
    edges         = hb.edges()         # np.array
    vmin          = hb.vmin()          # vtype
    vmax          = hb.vmax()          # vtype
    vtype         = hb.vtype()         # np.dtype
    equalbins     = hb.equalbins()     # bool

    limits        = hb.limits()        # np.array (vmin, vmax)
    binedges      = hb.binedges()      # np.array of size nbins+1 
    binedgesleft  = hb.binedgesleft()  # np.array of size nbins
    binedgesright = hb.binedgesright() # np.array of size nbins
    bincenters    = hb.bincenters()    # np.array of size nbins
    binwidth      = hb.binwidth()      # np.array of size nbins or scalar for Equal bins
    halfbinw      = hb.halfbinw()      # np.array of size nbins or scalar for Equal bins
    strrange      = hb.strrange(fmt)   # str of formatted vmin, vmax, nbins ex: 1-6-5

    # Print methods
    hb.print_attrs_defined()
    hb.print_attrs()
    hb.print_attrs_and_methods()

@see :py:class:`pyimgalgos.HBins`

This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@version $Id$

Created on 2016-01-15

@author Mikhail S. Dubrovin
"""
#------------------------------
__version__ = "$Revision$"
#------------------------------

#import math
import numpy as np

#------------------------------

class HBins() :
    """Hystogram-style bin parameters holder
    """
    def __init__(self, edges, nbins=None, vtype=np.float32):
        """Class constructor for
           - equal bins,       ex: hb = HBins((1,6), nbins=5)
           - or variable bins, ex: hb = HBins((1,2,4,6,10))

           Parameters:
           - edges - sequence of two or more bin edges
           - nbins - (int) number of bins for equal size bins
           - vtype - numpy type of bin values (optional parameter)          
        """
        self._name       = self.__class__.__name__
        self._edges      = np.array(edges, dtype=vtype)
        self._vmin       = min(self._edges)
        self._vmax       = max(self._edges)
        self._nbins      = int(nbins if nbins is not None else len(edges)-1)
        self._vtype      = vtype
        self._equalbins  = len(edges)==2 and nbins is not None

        self._limits     = None
        self._binwidth   = None 
        self._halfbinw   = None 
        self._binedges   = None 
        self._bincenters = None 
        self._inds       = None 
        self._indedges   = None 
        self._indcenters = None 
        self._strrange   = None


    def edges(self) :
        """Returns input sequence of edges"""
        return self._edges


    def vmin(self) :
        """Returns minimal value of the range"""
        return self._vmin


    def vmax(self) :
        """Returns miximal value of the range"""
        return self._vmax


    def nbins(self) :
        """Returns number of bins"""
        return self._nbins


    def vtype(self) :
        """Returns npumpy datatype for bin values"""
        return self._vtype


    def equalbins(self) :
        return self._equalbins


    def limits(self) :
        """Returns np.array of two ordered limits (vmin, vmax)"""
        if self._limits is None :
            self._limits = np.array((self._vmin, self._vmax), dtype=self._vtype)
        return self._limits


    def binedges(self) :
        """Returns np.array of nbins+1 values of bin edges"""
        if self._binedges is None : 
            if self._equalbins :
                self._binedges = np.linspace(self._vmin, self._vmax, self._nbins+1, endpoint=True, dtype=self._vtype)
            else :
                self._binedges = self._edges
        return self._binedges


    def binedgesleft(self) :
        """Returns np.array of nbins values of bin left edges"""
        return self.binedges()[:-1]


    def binedgesright(self) :
        """Returns np.array of nbins values of bin right edges"""
        return self.binedges()[1:]


    def binwidth(self) :
        """Returns np.array of nbins values of bin widths"""
        if self._binwidth is None :
            if self._equalbins :
                self._binwidth = float(self._vmax-self._vmin)/self._nbins
            else :
                self._binwidth = self.binedgesright() - self.binedgesleft()
        return self._binwidth


    def halfbinw(self) :
        """Returns np.array of nbins values of bin half-widths"""
        if self._halfbinw is None :
                self._halfbinw = 0.5 * self.binwidth()
        return self._halfbinw


    def bincenters(self) :
        """Returns np.array of nbins values of bin centers"""
        if self._bincenters is None :
            self._bincenters = self.binedgesleft() + self.halfbinw()
        return self._bincenters


    def strrange(self, fmt='%.0f-%.0f-%d') :
        """Returns string of range parameters"""
        if self._strrange is None :
            self._strrange =fmt % (self._vmin, self._vmax, self._nbins)
        return self._strrange


    def print_attrs(self) :
        print 'Attributes of the %s object' % self._name
        for k,v in self.__dict__.items() :
            print '  %s : %s' % (k.ljust(16), str(v))


    def print_attrs_defined(self) :
        print 'Attributes (not None) of the %s object' % self._name
        for k,v in self.__dict__.items() :
            if v is None : continue
            print '  %s : %s' % (k.ljust(16), str(v))


    def print_attrs_and_methods(self) :
        print 'Methods & attributes of the %s object' % self._name
        for m in dir(self) :
            print '  %s' % (str(m).ljust(16))

#------------------------------

def test(o, cmt='') :

    print '%s\n%s\n' % (80*'_', cmt)

    o.print_attrs_and_methods()
    o.print_attrs_defined()
    print 'nbins = %d' %   o.nbins()
    print 'limits',        o.limits()
    print 'binedges',      o.binedges()
    print 'binedgesleft',  o.binedgesleft()
    print 'binedgesright', o.binedgesright()
    print 'bincenters',    o.bincenters()  
    print 'binwidth',      o.binwidth()  
    print 'halfbinw',      o.halfbinw()  
    print 'strrange',      o.strrange()  
    o.print_attrs_defined()
    print '%s' % (80*'_')
    
#------------------------------

if __name__ == "__main__" :

    o1 = HBins((1,6), 5);      test(o1, 'Test HBins for EQUAL BINS')
    o2 = HBins((1,2,4,6,10));  test(o2, 'Test HBins for VARIABLE BINS')

#------------------------------
