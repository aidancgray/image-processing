# process_darks.py
# 08/04/2023
# Aidan Gray
# aidan.gray@idg.jhu.edu
#
# This is a script to average a bunch of dark fits images

from astropy.io import fits
import os
import sys
import glob
import numpy as np


def get_fits_size(filename):
    fitsFile = fits.open(filename)
    fitsData = fitsFile[0].data
    dataSize = fitsData.shape

    return dataSize


def open_fits(filename):
    fitsFile = fits.open(filename)
    fitsData = fitsFile[0].data
    fitsHdr = fitsFile[0].header

    return fitsData, fitsHdr


def loop_thru_dir(filePath):
    directoryList = glob.glob(filePath+'raw-*')
    directoryList.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    print("Combining images in: "+filePath)

    numFiles = len(directoryList)

    if numFiles < 1:
        sys.exit("ERROR: Directory is empty")

    firstData, firstHdr = open_fits(directoryList[0])
    totalData = firstData
    combHdr = firstHdr

    if numFiles > 1:
        for filename in directoryList[1:]:
            print(f"Adding: {filename}")
            tmpData, tmpHdr = open_fits(filename)
            totalData = totalData + tmpData

    combData = np.int16(totalData / numFiles)

    return combData, combHdr


if __name__ == "__main__":
    filePath = sys.argv[1]
    if len(sys.argv) > 2:
        combFilename = sys.argv[2]
    else:
        combFilename = "combined-darks.fits"

    if filePath[0] == '~':
        filePath = os.path.expanduser('~')+filePath[1:]

    if filePath[len(filePath)-1] != '/':
        filePath = filePath+'/'

    if not os.path.exists(filePath):
        sys.exit("ERROR: That file path does not exist.")
    else:
        combData, combHdr = loop_thru_dir(filePath)
        fits.writeto(
            filename=filePath+combFilename, 
            data=combData, 
            header=combHdr,
            overwrite=True
            )
