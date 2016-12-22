# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 21:34:11 2016

Flattens all 48-bit RGB tiffs in a directory, and saves them in a new folder 
(with naming modifications to pad and include time in the name). 
Saved as 16-bit tifs

Functions:
flatten48bitimages : flattens all 16-bit RGB (48bit) images in a folder.

2016Dec21: Modified to try to be more system-independent by using os.path.join:
UPDATED VERSION NOT TESTED YET.

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
                        tifmod.imsave(os.path.join(mydirname, newfoldername,
                                      newfilename), currentimage,
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
    if os.path.exists(os.path.join(mydirname, mynewfolder)):
        for k in range(maxk):
            if not os.path.exists(os.path.join(mydirname, mynewfolder, 
                                               str(k))):
                print(mynewfolder + ' already exists in this directory.')
                response = input(
                    'Make new folder ' + os.path.join(mynewfolder, str(k)) +
                    '? y/n')
                if response is 'y' or response is 'Y':
                    mynewfolder = os.path.join(mynewfolder, str(k))
                    os.mkdir(os.path.join(mydirname, mynewfolder))
                    break
                elif response is 'n' or response is 'N':
                    mynewfolder = []
                    break
                else:
                    print('Invalid response, treated as n')
                    mynewfolder = []
                    break
            elif os.path.exists(os.path.join(mydirname, mynewfolder, str(k))
                                ) and k == maxk-1:
                mynewfolder = []
                print('Failed to create folder with unique name in '
                    + str(maxk) + ' tries')
                break
            else:
                continue
            break
    else:
        os.mkdir(os.path.join(mydirname, mynewfolder))

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

#flatten48bitimages(zygdir)