#!/usr/bin/python3
# subtract_darks.py
# 6/12/2023
# Aidan Gray
# aidan.gray@idg.jhu.edu
#
# This is a script to process a batch of images

from astropy.io import fits
from matplotlib import pyplot as plt
import os
import sys
import numpy as np 
import glob
import PyGuide
import csv

# Switches ########################################
PIXEL_OUTPUT = True
DISPLAY_TARGETS = True
SUBTRACT_BIAS = False
SUBTRACT_DARK = False
BIAS_FILE = ''  # .fits filetype
DARK_FILE = ''  # .fits filetype
SAVE_PROCESSED_IMAGE = False
######################################################

# Constants #######################################
ZERO_PIXEL_DELTA = [0, 0]
ZERO_PIXEL = [1375, 1100]  # center of 2750x2200
PIXEL_SIZE = 0.00454  # mm
MASK_X_LEFT = 0  # 1300
MASK_X_RIGHT = 2749  # 1450
MASK_Y_BOTTOM = 990  # 1050
MASK_Y_TOP = 1190  # 1150
######################################################

# CCD Parameters for PyGuide init #################
BIAS_LEVEL = 0  # subtraction done using bias image
GAIN = 0.27  # e-/ADU
READ_NOISE = 3.5  # e-
MAX_COUNTS = 65536
FULL_WELL = 17000
######################################################


def write_to_csv(dataFile, dataList):
    print("Writing data to "+dataFile)
    with open(dataFile, 'w', newline='') as dF:
        wr = csv.writer(dF, dialect='excel', delimiter=',')
        wr.writerow([
            'filename', 
            'x (px)', 
            'y (px)', 
            'expTime (s)', 
            'flux', 
            'counts', 
            'fwhm (px)',
            ])
        for imageData in dataList:
            wr.writerows(imageData)
    print("Done")


def pyguide_checking(imgArray):
    """
    Uses PyGuide to find stars, get counts, and determine if new exposure time
    is necessary.

    Input:
    - imgArray  numpy array from the CCD

    Output:
    - True if exposure was good, False if bad
    - True if exposure time should be decreased, False if increased
    """
    # search image for stars
    centroidData, imageStats = PyGuide.findStars(
        imgArray,
        mask=maskArray,
        satMask=None,
        ccdInfo=CCDInfo,
        )

    # keep track of targets
    goodTargets = []
    lowTargets = 0
    highTargets = 0

    print(f"these are the {len(centroidData)} stars pyguide found:")
    bestTarget = None
    for centroid in centroidData:
        # for each star, measure its shape
        shapeData = PyGuide.starShape(
            np.asarray(imgArray, dtype="float32"),  # had to explicitly cast
            mask=None,
            xyCtr=centroid.xyCtr,
            rad=centroid.rad,
            )
        if not shapeData.isOK:
            print("starShape failed: %s" % (shapeData.msgStr,))
        else:
            print("xyCenter=[%.2f, %.2f] CCD Pixel Counts=%.1f, FWHM=%.1f, "
                  "BKGND=%.1f, chiSq=%.2f" %
                  (centroid.xyCtr[0], centroid.xyCtr[1], shapeData.ampl,
                   shapeData.fwhm, shapeData.bkgnd, shapeData.chiSq))
            # if shapeData.ampl < 0.2*MAX_COUNTS:
            #    lowTargets+=1
            # elif shapeData.ampl > 0.9*MAX_COUNTS:
            #    highTargets+=1
            # else:
            #    goodTargets.append([centroid,shapeData])
            if bestTarget is None:
                bestTarget = [centroid, shapeData]
            else:
                if shapeData.ampl >= bestTarget[1].ampl:
                    bestTarget = bestTarget = [centroid, shapeData]
            goodTargets.append([centroid, shapeData])
    print()

    print(str(len(goodTargets))+" targets are in the linear (20-90%) range ---"
          + " "+str(lowTargets)+" low targets --- "+str(highTargets) + 
          " high targets")

    # highlight detections
    # size of green circle scales with total counts
    # bigger circles for brigher stars
    if DISPLAY_TARGETS:
        plt.clf()
        # vmin/vmax help with contrast
        plt.imshow(imgArray, cmap="gray", vmin=200, vmax=MAX_COUNTS)
        plt.ion()
        plt.show()
        for centroid in centroidData:
            # offset by half a pixel to match imshow with 0,0 at pixel center 
            # rather than edge
            xyCtr = centroid.xyCtr + np.array([-0.5, -0.5])
            counts = centroid.counts
            plt.scatter(xyCtr[0], xyCtr[1], s=counts/MAX_COUNTS, 
                        marker="o", edgecolors="lime", facecolors="none",)
        plt.gca().invert_yaxis()
        plt.draw()
        plt.pause(0.1)

    # Successful exposure, return True. The False is thrown away
    # if len(goodTargets) >= 1:
    #     goodTargets = [goodTargets[0]]
    # if bestTarget != None:
    #     return [bestTarget]
    # else:
    #     return goodTargets
    return goodTargets


def single_image(fileName):
    """
    Function to process a single raw FITS.

    Input:
    - fileName      Name of absolute path to the raw FITS file

    Output:
    - dataList      List of coordinate points & data
    """
    rawFile = fits.open(fileName)
    rawData = rawFile[0].data
    rawHdr = rawFile[0].header

    dataList = []
    expTime = rawHdr['EXPTIME']
    
    if SUBTRACT_BIAS:
        biasFile = fits.open(BIAS_FILE)
        biasData = biasFile[0].data
        prcData = np.subtract(rawData, biasData)
    else:
        prcData = rawData

    goodTargets = pyguide_checking(prcData)

    if len(goodTargets) > 0:
        # dataList.append([fileName])
        for target in goodTargets:
            xPixel = target[0].xyCtr[0]
            yPixel = target[0].xyCtr[1]

            fluxTarg = target[0].counts
            countsTarg = target[1].ampl
            fwhmTarg = target[1].fwhm
            # bkgndTarg = target[1].bkgnd
            # chiSqTarg = target[1].chiSq

            targetData = [
                fileName, 
                xPixel, 
                yPixel, 
                expTime, 
                fluxTarg, 
                countsTarg, 
                fwhmTarg,
                ]
            dataList.append(targetData)

    if SAVE_PROCESSED_IMAGE:
        tmpFileName = fileName.split('raw')[1]
        prcFileName = 'prc'+tmpFileName
        fits.writeto(filePath+prcFileName, prcData, rawHdr)

    return dataList


def loop_thru_dir(filePath):
    """
    Function to loop through given directory and open all raw FITS.

    Input:
    - filePath      Name of the directory containing the raw FITS files

    Output:
    - dataList      List of coordinate points & data
    """
    dataList = []
    directoryList = glob.glob(filePath+'raw-*')
    directoryList.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    print("Processing images in: "+filePath)
    # print("List of filenames: "+repr(directoryList))
    
    for fileName in directoryList:
        print(f"Processing: {fileName}")
        dataListTemp = single_image(fileName)
        dataList.append(dataListTemp)

    return dataList


if __name__ == "__main__":
    filePath = sys.argv[1]
    dataFile = sys.argv[2]
    # dataFile = filePath[:-17]+filePath[-8:-4]+'csv'

    # filePath = input("Enter path to file (eg. ~/Pictures/SX_CCD): ")
    # dataFile = input("Enter desired output file (eg. data.csv): ")

    if filePath[0] == '~':
        filePath = os.path.expanduser('~')+filePath[1:]

    if dataFile[0] == '~':
        dataFile = os.path.expanduser('~')+dataFile[1:]

    CCDInfo = PyGuide.CCDInfo(
        bias=BIAS_LEVEL,  # image bias, in ADU
        readNoise=READ_NOISE,  # read noise, in e-
        ccdGain=GAIN,  # inverse ccd gain, in e-/ADU
        )

    maskArray = np.ones((2200, 2750))

    for row in range(MASK_Y_BOTTOM, MASK_Y_TOP+1):
        for col in range(MASK_X_LEFT, MASK_X_RIGHT+1):
            maskArray[row][col] = 0

    if filePath[len(filePath)-5:] == '.fits':
        dataListTemp = single_image(filePath)
        dataList = []
        dataList.append(dataListTemp)
        write_to_csv(dataFile, dataList)

    else:
        if filePath[len(filePath)-1] != '/':
            filePath = filePath+'/'

        if not os.path.exists(filePath):
            print("ERROR: That file path does not exist.")
            sys.exit()
        else:
            dataList = loop_thru_dir(filePath)
            write_to_csv(dataFile, dataList)
    if DISPLAY_TARGETS:
        input("Press ENTER to exit")
