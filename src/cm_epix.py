from __future__ import division
import numpy as np


def _neighborImg(img):
    img_up = np.roll(img, 1, axis=0)
    img_up[0, :] = 0
    
    img_down = np.roll(img, -1, axis=0)
    img_down[-1, :] = 0

    img_left = np.roll(img, 1, axis=1)
    img_left[:, 0] = 0

    img_right = np.roll(img, -1, axis=1)
    img_right[:, -1] = 0

    return np.amax(np.array([img_up, img_down, img_left, img_right]), axis=0)


def cm_epix(img, rms, maxCorr=30, histoRange=30, colrow=3,
            minFrac=0.25, normAll=False):
    '''
    Parameters:
    img - image on which the common mode is applied
    rms - array of noise values for each pixel with same shape as img
    maxCorr - (optional) maximum correction applied. If common mode
              correction is larger than this value, no correction will
              be applied
    histoRange - (optional) all pixels above this parameter are masked
    colrow - (optional) decides what is corrected. If 1, only the columns
             are corrected. If 2, only the rows are corrected. And if 3,
             both are corrected
    minFrac - (optional) the minimum fraction of pixels required to be
              left in a row or column after applying the mask and rejecting
              high pixels and their neighbors
    normAll - (optional) if true, will subtract the mean from the full
              image with the masked applied
    '''
    img = img.reshape(704, 768)
    imgThres = np.copy(img)
    imgThres[img >= rms * 10] = 1
    imgThres[img < rms * 10] = 0
    imgThres += _neighborImg(imgThres)
    imgThres += imgThres + (abs(img) > histoRange)

    maskedImg = np.ma.masked_array(img, imgThres)
    if normAll:
      maskedImg -= maskedImg.mean()

    if colrow % 2 == 1:
        rs = maskedImg.reshape(704 // 2 , 768 * 2 , order='F')
        rscount = np.ma.count_masked(rs, axis=0)
        rsmed = np.ma.median(rs, axis=0)
        rsmed[abs(rsmed) > maxCorr] = 0
        rsmed[rscount > ((1. - minFrac) * 352)] = 0
        imgCorr = np.ma.masked_array((rs.data -
                    rsmed[None, :]).data.reshape(704, 768, order='F'), imgThres)
    else:
        imgCorr = maskedImg.copy()

    if colrow >= 2:
        rs = imgCorr.reshape(704 * 8, 96)
        rscount = np.ma.count_masked(rs, axis=1)
        rsmed = np.ma.median(rs, axis=1)
        rsmed[abs(rsmed) > maxCorr]=0
        rsmed[rscount > ((1. - minFrac) * 96)] = 0
        imgCorr = (np.ma.masked_array(rs.data -
                        rsmed[:, None], imgThres)).reshape(704, 768)

    img[:] = imgCorr.data
