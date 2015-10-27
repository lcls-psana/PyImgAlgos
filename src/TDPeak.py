#!/usr/bin/env python

#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  class TDPeak
#
#------------------------------------------------------------------------

"""TDPeak - this class inherits from PeakData information holder/accessor class.

Works together with TDFileContainer and TDGroup classes.


This software was developed for the LCLS project.
If you use all or part of it, please give an appropriate acknowledgment.

@see TDFileContainer - loads/holds text data from class and provides per-event-indexed access. 
@see TDGroup - holds a list of records associated with a group of records.

@version $Id$

@author Mikhail S. Dubrovin
"""
#------------------------------
__version__ = "$Revision$"
# $Source$
##-----------------------------

from pyimgalgos.PeakData import PeakData

##-----------------------------

class TDPeak(PeakData) :
    """User-defined object for text data record parsing/processind/access.
    """
    def __init__(self, rec) :
        """Constructor
        @param rec - text data record from file
        """
        PeakData.__init__(self, rec)

    # Interface methods used in TDEvent
    def print_short(self) :
        self.print_peak_data_short()

##-----------------------------
##-----------------------------
##-----------------------------
