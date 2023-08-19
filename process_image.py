#!/usr/bin/python3
# geb_img_prc.py
# 8/18/2023
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
import shlex
import argparse

# Constants #######################################
ZERO_PIXEL_DELTA = [0, 0]
ZERO_PIXEL = [1375, 1100]  # center of 2750x2200
PIXEL_SIZE = 0.00454  # mm
MASK_X_LEFT = 0 
MASK_X_RIGHT = 2749 
MASK_Y_BOTTOM = 0 
MASK_Y_TOP = 2199  
######################################################

# CCD Parameters for PyGuide init #################
BIAS_LEVEL = 0  # subtraction done using bias image
GAIN = 0.27  # e-/ADU
READ_NOISE = 3.5  # e-
MAX_COUNTS = 65536
FULL_WELL = 17000
######################################################

CCD_INFO = PyGuide.CCDInfo(
        bias=BIAS_LEVEL,  # image bias, in ADU
        readNoise=READ_NOISE,  # read noise, in e-
        ccdGain=GAIN,  # inverse ccd gain, in e-/ADU
        )

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def open_fits(fileName):
    fitsFile = fits.open(fileName)
    fitsData = fitsFile[0].data
    fitsHdr = fitsFile[0].header

    return fitsData, fitsHdr


def write_fits(filepath, filename, data, hdr):
    fits.writeto(filepath+filename, data, hdr, overwrite=True)


def single_image_crop(filename, data, hdr, psSize, ctrs):
    xL_t = int(round(ctrs['t']['x']) - (psSize/2))
    yL_t = int(round(ctrs['t']['y']) - (psSize/2))
    xH_t = xL_t + psSize 
    yH_t = yL_t + psSize
    
    xL_c = int(round(ctrs['c']['x']) - (psSize/2))
    yL_c = int(round(ctrs['c']['y']) - (psSize/2))
    xH_c = xL_c + psSize 
    yH_c = yL_c + psSize
    
    xL_b = int(round(ctrs['b']['x']) - (psSize/2))
    yL_b = int(round(ctrs['b']['y']) - (psSize/2))
    xH_b = xL_b + psSize 
    yH_b = yL_b + psSize
    
    prcData_t = data[yL_t:yH_t, xL_t:xH_t]
    prcData_c = data[yL_c:yH_c, xL_c:xH_c]
    prcData_b = data[yL_b:yH_b, xL_b:xH_b]
    
    tmpFilename = filename.split('dsi-')[1]
    tmpFilepath = filename.split('dsi-')[0]

    prcFilename_t = 'post/t-'+tmpFilename
    prcFilename_c = 'post/c-'+tmpFilename
    prcFilename_b = 'post/b-'+tmpFilename
    
    write_fits(tmpFilepath, prcFilename_t, prcData_t, hdr)
    write_fits(tmpFilepath, prcFilename_c, prcData_c, hdr)
    write_fits(tmpFilepath, prcFilename_b, prcData_b, hdr)


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
    return goodTargets[:3]


def process(fileName, darkFileName, numTargets, fwhmBool, rmsBool, psSize):
    rawData, rawHdr = open_fits(fileName)

    dataList = []
    expTime = rawHdr['EXPTIME']
    
    if darkFileName is not None:
        darkData, darkHdr = open_fits(darkFileName)
        prcData = np.subtract(rawData, darkData)

    else:
        prcData = rawData

    goodTargets = pyguide_checking(prcData)

    if len(goodTargets) > 0:
        for target in goodTargets:
            xPixel = target[0].xyCtr[0]
            yPixel = target[0].xyCtr[1]

            fluxTarg = target[0].counts
            countsTarg = target[1].ampl

            if fwhmBool:
                fwhmTarg = target[1].fwhm
            else:
                fwhmTarg = '-'

            targetData = [fileName, 
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
        fits.writeto(filePath+prcFileName, prcData, rawHdr, overwrite=True)

    if SAVE_POSTAGE_IMAGES:
        tmpCtrs = {'t': {'x': -1, 'y': -1}, 
                   'c': {'x': -1, 'y': -1}, 
                   'b': {'x': -1, 'y': -1},
                   }

        for spot in dataList:
            if spot[2] > 1000 and spot[2] < 1200:
                if spot[1] < 400:
                    tmpCtrs['t']['x'] = spot[1]
                    tmpCtrs['t']['y'] = spot[2]
                elif spot[1] > 2350:
                    tmpCtrs['b']['x'] = spot[1]
                    tmpCtrs['b']['y'] = spot[2]
                elif spot[1] > 1300 and spot[1] < 1450:
                    tmpCtrs['c']['x'] = spot[1]
                    tmpCtrs['c']['y'] = spot[2]


        single_image_crop(filename=fileName, 
                          data=prcData,
                          hdr=rawHdr,
                          psSize=40,
                          ctrs=tmpCtrs,
                          )

    return dataList


def loop_thru_dir(filePath, wildcardSearch='raw-*'):
    """
    Function to loop through given directory and open all raw FITS.

    Input:
    - filePath      Name of the directory containing the raw FITS files

    Output:
    - dataList      List of coordinate points & data
    """
    dataList = []
    directoryList = glob.glob(filePath+wildcardSearch)
    directoryList.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    print("Processing images in: "+filePath)
    # print("List of filenames: "+repr(directoryList))
    
    for fileName in directoryList:
        print(f"Processing: {fileName}")
        dataListTemp = single_image(fileName)
        dataList.append(dataListTemp)

    return dataList

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if isinstance(argv, str):
        argv = shlex.split(argv)
    
    parser = argparse.ArgumentParser(sys.argv[0])
    parser.add_argument(
        '--filePath', 
        type=str, 
        default=None,
        help='path to the image file',
        )
    parser.add_argument(
        '--darkFile', 
        type=str, 
        default=None,
        help='(optional) path to the image file ' \
                'used for dark subtraction',
        )
    parser.add_argument(
        '--numTargets', 
        type=int, 
        default=-1,
        help='(optional) number of expected targets',
        )
    parser.add_argument(
        '--FWHM', 
        type=bool, 
        default=False,
        help='bool: measures target FWHM if True',
        )
    parser.add_argument(
        '--RMS', 
        type=bool, 
        default=False,
        help='bool: measures target RMS if True',
        )
    parser.add_argument(
        '--psSize', 
        type=int, 
        default=-1,
        help='(optional) save postage stamp images ' \
                'of the target with these square dimensions',
        )

    opts = parser.parse_args(argv)
    
    filePath = opts.filePath
    darkFile = opts.darkFile
    numTargets = opts.numTargets
    fwhmBool = opts.fwhm
    rmsBool = opts.RMS
    psSize = opts.psSize

    if filePath is None:
        sys.exit('ERROR: FITS file not specified')
    if filePath[0] == '~':
        filePath = os.path.expanduser('~') + filePath[1:]
    fileName = os.path.basename(os.path.normpath(filePath))
    
    if darkFile is not None:
        if darkFile[0] == '~':
            darkFile = os.path.expanduser('~') + darkFile[1:]

    maskArray = np.ones((2200, 2750))

    for row in range(MASK_Y_BOTTOM, MASK_Y_TOP+1):
        for col in range(MASK_X_LEFT, MASK_X_RIGHT+1):
            maskArray[row][col] = 0

    if filePath[len(filePath)-5:] == '.fits':
        dataFile, dataListTemp = process(
            filePath, darkFile, numTargets, fwhmBool, rmsBool, psSize)
        dataList = []
        dataList.append(dataListTemp)
        write_to_csv(dataFile, dataList)
    else:
        sys.exit('ERROR: FITS file not specified')

'''
    else:
        if filePath[len(filePath)-1] != '/':
            filePath = filePath+'/'

        if not os.path.exists(filePath):
            print("ERROR: That file path does not exist.")
            sys.exit()
        else:
            if len(sys.argv) > 3:
                dataList = loop_thru_dir(filePath, sys.argv[3])
            else:
                dataList = loop_thru_dir(filePath)
            write_to_csv(dataFile, dataList)
    if DISPLAY_TARGETS:
        input("Press ENTER to exit")
'''

if __name__ == "__main__":
    main()