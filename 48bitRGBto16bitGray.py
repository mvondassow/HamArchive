# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 21:34:11 2016

Flattens all 48-bit RGB tiffs in a directory, and saves them in a new folder 
(with naming modifications to pad and include time in the name). 
Saved as 16-bit tifs

Functions:
flatten48bitimages : flattens alls images in 

@author: Michelangelo
"""
import numpy as np
import skimage.external.tifffile as tifmod

import os
import time
import glob


def flatten48bitimages(mydirname, **kwargs):
    """
    Convert all 48-bit RGB tiffs to 16-bit grayscale tifs. Averages three color
    channels, and saves in file.
    
    params
    -------
    mydirname: string
        name of directory containing images
    kwargs: 
        myextension: string extension (currently only handles tif or tiff)

    returns
    ------
    nothing
    """
    # 'extension' kwarg as file extension
    if (kwargs.get('extension') is None):
        myextension = 'tiff'
    else:
        myextension = '.' + kwargs.get('extension')

    # Create generic pathname for images
    mypathname = mydirname + '*.' + myextension
    # Get list of file names matching pattern along path
    myfilelist = glob.glob(mypathname)
    # Check that some files exist that fit glob name; then create new folder.
    # and flatten and save images.
    if len(myfilelist)>0:
        # create a new folder in path and return its name
        newfoldername = mynewfolder(mydirname, 'flattened', maxk=10)
        if len(newfoldername)>0:
            # decide how much padding to add to names.
            paddinglength = len(str(len(myfilelist)))
            # Run through each image: create new name, flatten channels, and
            # save with new name.
            for k in range(len(myfilelist)):
                # Get image
                currentfile = myfilelist[k]
                currentimage = tifmod.imread(myfilelist[k])
                # get file modification date (SEE os module notes about
                # st_*time: st_ctime seems to give time when file was copied, 
                # not when first created).
                filedate = os.stat(currentfile).st_mtime
                # Make new file name, padding to max length of old file name.
                # Also adds date/time in serial date up to seconds
                newfilename = currentfile[len(mydirname):][:-(
                                        len(myextension)+1)].zfill(
                                        paddinglength) + "_" + str(
                                        round(filedate)) + ".tif"
                # check if image reads as numpy.ndarray and has correct type 
                # (16 bit) and shape (3 channels).
                if type(currentimage) is numpy.ndarray:
                    if (currentimage.dtype == 'uint16') and (
                            currentimage.ndim == 3):
                        # Flatten color channels to 16-bit
                        currentimage = np.mean(
                            currentimage, 2).astype(np.uint16)
                        # For the life of me, I can't get this module to save
                        # metadata as it claims it will. Just saving in
                        # description.
                        tifmod.imsave(mydirname + newfoldername + '\\' +
                                      newfilename,
                                      currentimage,
                                      description=currentfile)
                    else:
                        print(currentfile + ' is not of the right type.')
                        
                else:
                        print(currentfile + ' not read.')
        else:
            print('Folder not made, so new images not saved')
    else:
        print('No files with correct extension found.')


def mynewfolder(mydirname, mynewfolder, maxk = 10):
    """
    Create a new directory in path 'mydirname. with name 'mynewfolder'. If 
    folder with name 'mynewfolder' already exists in 'mydirname', it asks to
    create a new one, incrementing from 0 to maxk; stops if it reaches maxk.
    
    Parameters:
    ----------
    mydirname : string
        valid path name
    mynewfolder : string
        valid folder name
    maxk : int
    
    Returns:
    ---------
    new folder name after incrementing (empty list - [] - if cannot create
    folder.)
    """
    if os.path.exists(mydirname + mynewfolder):
        for k in range(maxk):
            if not os.path.exists(mydirname + mynewfolder + str(k)):
                print(mynewfolder + ' already exists in this directory.')
                response = input(
                    'Make new folder ' + mynewfolder + str(k) + '? y/n')
                if response is 'y' or response is 'Y':
                    mynewfolder = mynewfolder + str(k)
                    os.mkdir(mydirname + mynewfolder)
                    break
                elif response is 'n' or response is 'N':
                    mynewfolder = []
                    break
                else:
                    print('Invalid response, treated as n')
                    mynewfolder = []
                    break
            elif os.path.exists(
                    mydirname + mynewfolder + str(k)) and k == maxk-1:
                mynewfolder = []
                print('Failed to create folder with unique name in '
                    + str(maxk) + ' tries')
                break
            else:
                continue
            break
    else:
        os.mkdir(mydirname + mynewfolder)

    return mynewfolder

#zygdir = 'C:\\Users\\Michelangelo\\Documents\\Ham\\CellVolumeRegulation_Zygotes\\'
#myfolders = []
#pathnames = glob.glob(zygdir + '*')
#for k in range(len(pathnames)):
#    if not pathnames[k].endswith('.tiff'):
#        myfolders = myfolders + [pathnames[k]]
#        
#for k in range(len(myfolders)):
#    flatten48bitimages(myfolders[k]+'\\')

'''
# Steps to process files from ImageJ results
1) Get data in.
import csv
myfolder = 'C:\\Users\\Michelangelo\\Documents\\Ham\\CellVolumeRegulation_Zygotes\\CVR_rib06_2016Aug03_12_00_58AM\\flattened\\'
myfile = 'Results_CVR_rib06_MeasEmbV5.xls'
with open(myfolder+myfile) as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
       ... (e.g. print(row[' ']) or print(row['Label']))
 Rows as dicts with keys from header lines.
with open(myfolder+myfile) as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    for row in reader:
       ...
 Gives rows as lists, including header rows.

2) Combine the following steps into one loop.
 2a: Split 'Label' to get time stamps and file codes
 2b: Identify XY centers based on 'Feret' 'FeretX' 'FeretY' and 'Angle'
 2c: calculate volume
 2d: Divide data into two sets (before and after transition) and delete frames
 during transition (identify manually)

3) Match ROIs from one frame to next based on overlap in XY center. In at
 least one case, an embryo was missed in one frame (crossed edge) so join if
 overlap is close between every other frame as well. Maybe use a clustering
 function?

4) check that </=1 embryo per cluster per frame: else, stop and report error.
  and group data by embryo (cluster) to create numpy arrays of time and volume
  for each embryo

 http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.vq.kmeans2.html#scipy.cluster.vq.kmeans2
'''
'''
import csv
myfolder = 'C:\\Users\\Michelangelo\\Documents\\Ham\\CellVolumeRegulation_Zygotes\\CVR_rib06_2016Aug03_12_00_58AM\\flattened\\'
myfile = 'Results_CVR_rib06_MeasEmbV5.xls'
with open(myfolder+myfile) as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    datalist = list(reader)

for row in datalist:
    # Partition label into sub strings and convert time and frame portions
    # to integers.
    row['Label'], row['Time'] = row['Label'].rsplit('_', 1)
    row['Label'], row['Image'] = row['Label'].rsplit(':', 1)
    row['Time'] = int(row['Time'])
    row['Image'] = int(row['Image'])

'''

'''
Looks like pandas might have all the functionality I want built in.
'''
import pandas
import numpy as np
import scipy.cluster.vq as vq

# File info
myfolder = 'C:\\Users\\Michelangelo\\Documents\\Ham\\CellVolumeRegulation_Zygotes\\CVR_rib06_2016Aug03_12_00_58AM\\flattened\\'
myfile = 'Results_CVR_rib06_MeasEmbV5.xls'
MicronsPerPixel = 0.338774005
TransitionImage = {'end1': 2, 'begin2': 3} # QCam saves from 0 to n-1. 'end1' 
# gives last frame (in QCams numbering) of the first part of the data series
# (before transition to high salinity); begin2 gives first frame of the 2nd
# part (after transition to high salinity).

# import data from file 
curdata = pandas.read_csv(myfolder+myfile, delimiter='\t')

# Split out frame (image) and time information from Label column 
curdata['Image'] = pandas.to_numeric(
    curdata['Label'].str.rsplit(':',1).str[1].str.rsplit('_',1).str[0])
curdata['Time'] = pandas.to_numeric(curdata['Label'].str.rsplit('_',1).str[1])

# add columns assigning membership to first or second part of data series.
curdata['part1'] = curdata['Image'] <= TransitionImage['end1']
curdata['part2'] = curdata['Image'] >= TransitionImage['begin2']

# get x and y position of embryo center (approximately)
curdata['xcent'] = curdata['FeretX'] + np.cos(
    curdata['FeretAngle']*2*np.pi/360)*curdata['Feret']/2
# ImageJ is treating y as increasing downwards, but a positive angle pointing
# upwards, so SUBTRACT sin(angle)*length.
curdata['ycent'] = curdata['FeretY'] - np.sin(
    curdata['FeretAngle']*2*np.pi/360)*curdata['Feret']/2
    
# Calculate volume in cubic micrometers.
curdata['Volume'] = (4*np.pi/3)*(
    MicronsPerPixel*(curdata['Feret']+curdata['MinFeret'])/2)**3

#TURNS OUT THAT IMAGEJ DOES NOT CALCULATE FERET ANGLE STARTING FROM FeretX AND
# FeretY: Either redo measurements from ROIs, adding centroid, or need to 
# figure out how it how it decides. FeretX and FeretY. Seems to flip enough
# that messes up clustering