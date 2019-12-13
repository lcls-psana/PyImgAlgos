from __future__ import print_function
from __future__ import division
import numpy as np
from skbeam.core.accumulators.histogram import Histogram
from PSCalib.CalibFileFinder import make_calib_file_name, find_calib_file
from PSCalib.h5constants import save as h5save
from PSCalib.h5constants import load as h5load

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

calib_type = 'photon_position'


class CentroidSmootherCalib(object):
    '''
    A calibration object that needs only to run once. Data and constants will
    be saved to be used by the CentroidSmoother object.
    '''
    def __init__(self, nBins=100):
        '''
        Parameters: nBins - (Optional) The number of bins for the histogramming
                            of the absolute distance from the center of the
                            pixel for each centroid
        '''
        self.nBins = nBins
        self.absDist = Histogram((self.nBins, 0.0, 0.5))


    def add(self, centroids):
        '''
        Parameters: An array of floats of the rows and columns of peak
                    centroids in the shape (x, 2) for x peaks.

        Computes, for each centroid, the absolute distance from the center of
        the pixel that the centroid is in. The values range from 0 to 0.5.
        '''
        flattened = centroids.flatten()
        absPosInPixel = abs(flattened - np.rint(flattened))
        self.absDist.fill(absPosInPixel)


    def save(self, ds, det, rNumBegin, rNumEnd=None):
        '''
        Parameters: ds - datasource as obtained via psana.DataSource()
                    det - detector used as obtained via psana.Detector()
                    rNumBegin - run number that the saved data begins with
                    rNumEnd - (Optional) run number that the saved data
                              ends with
        
        From the given parameters, saves the data that has been added by
        add() along with the bins which are computed from the value for
        nBins.
        '''
        cdir = ds.env().calibDir()
        src = str(det.name)

        absDist = comm.reduce(self.absDist.values)

        if rank == 0:
            nSum = absDist.cumsum()
            nTotal = absDist.sum()

            self.startBin = np.append(0, nSum[:-1] / (2 * nTotal))
            self.endBin = nSum / (2 * nTotal)

            fout  = make_calib_file_name(cdir, src, calib_type, rNumBegin, rNumEnd)
            writedict = {'version' : 1,
                         'startBin' : self.startBin,
                         'endBin' : self.endBin}
            h5save(fout, writedict)


class CentroidSmoother(object):
    '''
    With the datasource, detector and run number, this object will fetch
    the relevant data as saved by the CentroidSmootherCalib object and
    read it in.
    '''
    def __init__(self, ds, det, rNum):
        '''
        Parameters: ds - datasource as obtained via psana.DataSource()
                    det - detector used as obtained via psana.Detector()
                    rNum - run number of data
        '''
        cdir = ds.env().calibDir()
        src = str(det.name)
        fin = find_calib_file(cdir, src, calib_type, rNum)
        readdict = h5load(fin)

        self.startBin = readdict['startBin']
        self.endBin = readdict['endBin']
        self.nBins = len(self.endBin)


    def getSmoothedCentroids(self, offsets):
        '''
        Parameters: An array of floats the rows and columns of peak centroids
                    in the shape (x, 2) for x peaks.
        Output: An array of floats of the rows and columns of the smoothed
                centroids in the same shape as the input.

        Given an array of centroids, this function smooths and returns them.
        '''
        fracOffsets = offsets - np.rint(offsets)
        fracROffsets, fracCOffsets = fracOffsets[:, 0], fracOffsets[:, 1]

        rCentroid = self._correctedOffset(fracROffsets)
        cCentroid = self._correctedOffset(fracCOffsets)

        return np.dstack((rCentroid, cCentroid))[0]


    def _correctedOffset(self, offsets):
        signs = np.sign(offsets)
        signs[signs == 0] = np.random.choice([-1, 1], np.sum(signs == 0))
        offsets[abs(offsets) >= 0.5] = signs[abs(offsets) >= 0.5] * (0.5 - 1e-6)

        bins = (abs(offsets) / 0.5 * self.nBins).astype(int)
        randsInBin = np.random.uniform(0, 1, len(offsets)) * \
                     (self.endBin[bins] - self.startBin[bins])

        centroids = signs * (self.startBin[bins] + randsInBin)

        return centroids


if __name__ == "__main__":
    from ImgAlgos.PyAlgos import PyAlgos
    from psana import *
    import matplotlib.pyplot as plt

    ds = DataSource('exp=xcs06016:run=37:smd')
    det = Detector('epix_2')

    setOption('calib-dir', '/reg/neh/home/jscott/calib')

    alg = PyAlgos()
    alg.set_peak_selection_pars(npix_min=1, npix_max=4, amax_thr=0,
                                atot_thr=35, son_min=0)

    centroids = []

    for nevent, evt in enumerate(ds.events()):
        if nevent == 25:
            break

        nda = det.calib(evt)

        peaks = alg.peak_finder_v1(nda, thr_low=10, thr_high=30, radius=1, dr=0)
        for p in peaks:
            centroids.append([p[6], p[7]])

    print('%d total peaks.' % len(centroids))
    centroids = np.array(centroids)

    calib = CentroidSmootherCalib()
    calib.add(centroids)
    calib.save(ds, det, 0)

    cs = CentroidSmoother(ds, det, evt.run())

    smoothCents = cs.getSmoothedCentroids(centroids)


    smoothDist, edges = np.histogram(smoothCents.flatten(),
                                     bins=50, range=(-0.5, 0.5))
    smoothRD, edgesR = np.histogram(smoothCents[:, 0],
                                    bins=50, range=(-0.5, 0.5))
    smoothCD, edgesC = np.histogram(smoothCents[:, 1],
                                    bins=50, range=(-0.5, 0.5))

    plt.figure()
    plt.plot(calib.absDist.values)

    plt.figure()
    plt.plot(edgesR[:-1], smoothRD, marker='x')
    plt.title('Rows')
    plt.xlim(-0.5, 0.5)

    plt.figure()
    plt.plot(edgesC[:-1], smoothCD, marker='x')
    plt.title('Columns')
    plt.xlim(-0.5, 0.5)

    plt.figure()
    plt.plot(edges[:-1], smoothDist, marker='*')
    plt.xlim(-0.5, 0.5)

    plt.show()
