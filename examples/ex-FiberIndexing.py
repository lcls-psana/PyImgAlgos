#!/usr/bin/env python
#------------------------------

"""Uses global methods from pyimgalgos.FiberIndexing to
   - generate lattice in 3-d and reciprocal space
   - evaluate lattice node parameters sorted by node radius
   - rotate crystal lattice in beta and omega, find nodes close to Evald's sphere and generate look-up table
"""
from __future__ import print_function

#------------------------------
from pyimgalgos.FiberIndexing import *
#------------------------------

def make_index_table() :

    from pyimgalgos.GlobalUtils import str_tstamp

    fname = 'lut-cxif5315-r0169-%s.txt' % (str_tstamp())
    fout = open(fname,'w')
    fout.write('# file name: %s\n' % fname)

    # Lattice parameters
    a, b, c = 18.36, 26.65, 4.81        # Angstrom
    alpha, beta, gamma = 90, 90, 77.17  # 180 - 102.83 degree
    hmax, kmax, lmax = 5, 5, 0          # size of lattice to consider

    a1, a2, a3 = triclinic_primitive_vectors(a, b, c, alpha, beta, gamma)
    b1, b2, b3 = reciprocal_from_bravias(a1, a2, a3)

    msg1 = '\n# Triclinic crystal cell parameters:'\
         + '\n#   a = %.2f A\n#   b = %.2f A\n#   c = %.2f A' % (a, b, c)\
         + '\n#   alpha = %.2f deg\n#   beta  = %.2f deg\n#   gamma = %.2f deg' % (alpha, beta, gamma)

    msg2 = '\n# 3-d space primitive vectors:\n#   a1 = %s\n#   a2 = %s\n#   a3 = %s' %\
           (str(a1), str(a2), str(a3))

    msg3 = '\n# reciprocal space primitive vectors:\n#   b1 = %s\n#   b2 = %s\n#   b3 = %s' %\
           (str(b1), str(b2), str(b3))

    rec = '%s\n%s\n%s\n' % (msg1, msg2, msg3)
    print(rec)
    fout.write(rec)

    #for line in triclinic_primitive_vectors.__doc__.split('\n') : fout.write('\n# %s' % line)

    # Photon energy
    Egamma_eV = 6003.1936             # eV SIOC:SYS0:ML00:AO541
    wavelen_nm = wavelength_nm_from_energy_ev(Egamma_eV)
    evald_rad = wave_vector_value(Egamma_eV)
    sigma_q = 0.001 * evald_rad

    rec  = '\n# photon energy = %.4f eV' % (Egamma_eV)\
         + '\n# wavelength = %.4f A' % (wavelen_nm*10)\
         + '\n# wave number/Evald radius k = 1/lambda = %.6f 1/A' % (evald_rad)\
         + '\n# sigma_q   = %.6f 1/A (approximately pixel size/sample-to-detector distance = 100um/100mm)' % (sigma_q)\
         + '\n# 3*sigma_q = %.6f 1/A' % (3*sigma_q)\
         + '\n# %s\n\n' % (89*'_')

    print(rec)
    fout.write(rec)

    #test_lattice()
    test_lattice(b1, b2, b3, hmax, kmax, lmax, cdtype=np.float32)

    lattice_node_radius(b1, b2, b3, hmax, kmax, lmax, cdtype=np.float32)
    lattice_node_radius(b1, b2, b3, hmax, kmax, 1, cdtype=np.float32)
    #return
    #------------------------------

    import matplotlib.pyplot as plt
    import pyimgalgos.GlobalGraphics as gg

    # bin parameters for q in units of k = Evald's sphere radius [1/A]
    bpq = BinPars((-0.2, 0.2), 800, vtype=np.float32)

    # bin parameters for omega [degree] - fiber rotation angle around axis
    bpomega = BinPars((0., 180.), 360, vtype=np.float32)
    
    # bin parameters for beta [degree] - fiber axis tilt angle
    #bpbeta = BinPars((10., 20.), 3, vtype=np.float32, endpoint=True)
    bpbeta = BinPars((15., 15.), 1, vtype=np.float32, endpoint=False)
 
    lut = make_lookup_table(b1, b2, b3, hmax, kmax, lmax, np.float32, evald_rad, sigma_q, fout, bpq, bpomega, bpbeta)

    fout.close()
    print('\nFile with lookup table is saved: %s' % fname)

    img_range = (bpq.vmin, bpq.vmax, bpomega.vmax, bpomega.vmin) 
    axim = gg.plotImageLarge(lut, img_range=img_range, amp_range=None, figsize=(15,13),\
                      title='Image', origin='upper', window=(0.05,  0.06, 0.94, 0.94))
    axim.set_xlabel('$q_{H}$ ($1/\AA$)', fontsize=18)
    axim.set_ylabel('$\omega$ (degree)', fontsize=18)

    arrhi = np.sum(lut,0)
    fighi, axhi, hi = gg.hist1d(bpq.binedges, bins=bpq.nbins+1, amp_range=(bpq.vmin, bpq.vmax), weights=arrhi,\
                                color='b', show_stat=True, log=False,\
                                figsize=(15,5), axwin=(0.05, 0.12, 0.85, 0.80),\
                                title=None, xlabel='$q_{H}$ ($1/\AA$)', ylabel='Intensity', titwin=None)
    gg.show()

    #import matplotlib.pyplot  as plt
    #plt.imshow(lut)
    #plt.show()
    
#------------------------------

if __name__ == "__main__" :
    make_index_table()

#------------------------------
