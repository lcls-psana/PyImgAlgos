
### #!/usr/bin/env python

from pyimgalgos.NDArrGenerators import *


def plot_image(img, img_range=None, amp_range=None, figsize=(12,10)) :
    import pyimgalgos.GlobalGraphics as gg
    axim = gg.plotImageLarge(img, img_range, amp_range, figsize)
    gg.show()


def ex_random_standard(ntest) :
    img = random_standard(shape=(40,60), mu=200, sigma=25)
    plot_image(img)


def ex_random_exponential(ntest) :
    img = random_exponential(shape=(40,60), a0=100)
    plot_image(img)


def ex_random_one(ntest) :
    img = random_one(shape=(40,60), dtype=np.float32)
    plot_image(img)


def ex_random_256(ntest) :
    img = random_256(shape=(40,60), dtype=np.uint8)
    plot_image(img)


def ex_random_xffffffff(ntest):
    img = random_xffffffff(shape=(40,60), dtype=np.uint32, add=0xff000000)
    plot_image(img)


def ex_aranged_array(ntest):
    img = aranged_array(shape=(40,60), dtype=np.uint32)
    plot_image(img)


def ex_add_ring(ntest):
    img = random_standard(shape=(500, 500), mu=0, sigma=10)
    add_ring(img, amp=50, row=230, col=230, rad=150, sigma=2)
    add_ring(img, amp=99, row=280, col=280, rad=200, sigma=3)
    plot_image(img)


def ex_add_random_peaks(ntest):
    img = random_standard(shape=(500, 500), mu=0, sigma=10)
    peaks = add_random_peaks(img, npeaks=10, amean=100, arms=50, wmean=1.5, wrms=0.3)
    print('peaks:')
    for i, (r0, c0, a0, sigma) in enumerate(peaks) :
        print('  %04d  row=%6.1f  col=%6.1f  amp=%6.1f  sigma=%6.3f' % (i, r0, c0, a0, sigma))
    plot_image(img)


if __name__ == "__main__" :
    import sys; global sys
    tname = sys.argv[1] if len(sys.argv) > 1 else '0'
    print(50*'_', '\nTest %s:' % tname)
    if   tname == '0': ex_random_standard(tname); ex_add_random_peaks(tname); ex_add_ring(tname)
    elif tname == '1': ex_random_standard(tname)
    elif tname == '2': ex_random_exponential(tname)
    elif tname == '3': ex_random_one(tname)
    elif tname == '4': ex_random_256(tname)
    elif tname == '5': ex_random_xffffffff(tname)
    elif tname == '6': ex_aranged_array(tname)
    elif tname == '7': ex_add_ring(tname)
    elif tname == '8': ex_add_random_peaks(tname)
    else: print('Not-recognized test name: %s' % tname)
    sys.exit('End of test %s' % tname)

# EOF
