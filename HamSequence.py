# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 19:56:12 2016

@author: Michelangelo

Group and analyze blobs (embryos) from ImageJ macro output, given information
in info file.

Requires:
•Tab delimited file (with appropriate format and labels) containing info needed
to identify breaks in sequence, or measurements to treat specially.
•Tab delimited file with appropriate columns and labels containing data for
each sequence
•Scale info (CURRENTLY CODED IN SCRIPT)
•Directory path for files (CURRENTLY CODED IN SCRIPT)

'''
# Steps to process files from ImageJ results
1) Get data in.

2) Combine the following steps into one loop.
 2a: Split 'Label' to get time stamps and file codes
 2b: calculate volume
 2c: Divide data into two sets (before and after transition) and delete frames
 during transition (identify manually)

3) Match ROIs from one frame to next based on overlap in XY center. In at
 least one case, an embryo was missed in one frame (crossed edge) so join if
 overlap is close between every other frame as well.
'''

"""
import pandas
import numpy as np
import json
import TrackPoints
import matplotlib.pyplot as plt

# Parent directory
parentdir = 'C:\\Users\\Michelangelo\\Documents\\Ham\\'
# Name of file with info about sequences
infofile = parentdir + 'CellVolumeRegulation.txt'

# Scale in images
MicronsPerPixel = 0.338774005


def seqprocess(infodf, ind, folder, scale):
    """
    Import and process data from ImageJ macro.
    Note: image sequences may be broken up into parts (e.g. if capture stalled
    and was restarted), but are all assumed to be broken into two portions
    (before and after media change) unless 'end1' set to -1.

    File IO :
    --------
    Reads a file at location specified by folder and infodf.
    This file must contain a tab-delimited text file with the following
    columns: ' ' (int. ImageJ measurement indices), 'Label' (str. in format
    [!:_]*:[!:_]*:[0-9]*_[0-9]* . The 3rd glob is parsed to the 'Image' column,
    the 4th is parsed to the 'Time' column), 'Area' (float), 'X' (float), 'Y'
    (float), 'Major' (float), 'Minor' (float), 'Angle' (float), 'Feret'
    (float), 'Slice' (int. slice in image stack), 'FeretX' (float), 'FeretY'
    (float), 'FeretAngle' (float), 'MinFeret' (float).

    Parameters :
    -----------
    infodf : Pandas data frame with info about the group of image sequences
        columns: 'Sequence' (str. name of sequence), 'Part' (int. part of
        sequence), 'SetDir' (str. directory for set/subset of sequences),
        'SequenceDir' (str. directory within SetDir for specific part of
        sequence), 'MyFile' (str. file name for data for each sequence part),
        'End1' (int. image name at end of first portion of sequence), 'Begin2'
        (image name at beginning of second portion of sequence),
        'UseButMoving' (str. in JSON format for list of frames that may require
        special tratment), 'TrackMethod' (str. what method to use to track
        blobs), 'MeasIndexJoins ' (str in JSON format for sets of measurements
        that should be manually joined, set to -1 if images with usable data
        where blobs cannot be joined), 'DeleteMeasInds' (str. in JSON format,
        measurements to delete), 'Notes' (str.)
    ind : row index to get values from in infodf
    folder : parent directory path for SetDir.
    scale : scale for calculating volume

    Returns :
    Pandas data frame. Columns are the same as the text file, but the ' ' is
    renamed to 'IJind', 'Label' is parsed to give 'Image' and 'Time' columns,
    and 'Volume' is calculated. The data frames from multi-part sequences may
    be combined depending on user input.
    """
    # Sequence directory name. Assumes data is in sub folder 'flattened'!
    seqfolder = folder + infodf.SetDir[ind
                    ] + '\\' + infodf.SequenceDir[ind] + '\\flattened\\'
    # Data file name
    seqfile = infodf.MyFile[ind]
    # QCam saves from 0 to n-1. 'end1' gives last frame (in QCams numbering) of
    # the first part of the data series (before transition to high salinity);
    # begin2 gives first frame of the 2nd part (after transition to high
    # salinity).
    end1 = int(infodf.End1[ind])
    begin2 = int(infodf.Begin2[ind])

    # Get list of measurements to remove, if any.
    delmeaslist = json.loads(infodf.DeleteMeasInds[ind])
    # Get list of images in which embryos are moving enough to upset tracking.
    movelist = json.loads(infodf.UseButMoving[ind])
    """
    For now, see if auto-tracking is sufficient for those sequences with moving
    blobs/embryos that can be tracked (e.g. marked 'Manual' or 'Auto', not
    'Cannot'). If so, don't need to parse metadata about moving embryos, etc.
    """

    # import data from file
    curdata = pandas.read_csv(seqfolder+seqfile, delimiter='\t')

    # Replace ImageJ index column name (ImageJ labels as ' ') with 'IJind'
    newcolnames = dict(zip(curdata.columns, curdata.columns))
    newcolnames[' '] = 'IJind'
    curdata.rename(columns=newcolnames, inplace=True)

    # Split out frame (image) and time information from Label column
    curdata['Image'] = pandas.to_numeric(
        curdata['Label'].str.rsplit(':', 1).str[1].str.rsplit('_', 1).str[0])
    curdata['Time'] = pandas.to_numeric(
                                curdata['Label'].str.rsplit('_', 1).str[1])

    # Calculate volume in cubic micrometers.
    curdata['Volume'] = (4*np.pi/3)*(
        scale*(curdata['Feret']+curdata['MinFeret'])/2)**3

    # Remove unusable rows.
    # ImageJ spits out column of indices starting from 1; because the data file
    # may be edited, Pandas' indices might not be simply related to ImageJ's.
    # Therefore, create list of indices (dataframe) for which value in column
    # 'IJinds' matches any value in delmeaslist, or matches image names between
    # end1 and begin2 (to get rid of unmeasurable images... this code could
    # definitely be optimized.
    dropinds = []
    for k in curdata.index:
        if any(x == curdata['IJind'][k] for x in delmeaslist) | (
                                    begin2 > curdata['Image'][k] > end1):
                dropinds.append(k)

    # Define column to specify if images are of blobs in first medium (0) or
    # second medium (1)
    curdata['Media'] = curdata['Image'] >= begin2

    # Define column for images that are moving.
    moveimages = []
    for item in movelist:
        if len(item) == 2:
            moveimages = moveimages + (list(range(item[0], item[1]+1)))
        elif (len(item) != 2) & (len(item) != 0):
            print(len(item))
            raise SystemExit(
                'Error: "UseButMoving" has invalid sub-list length.')

    curdata['Moving'] = curdata['Image'].isin(moveimages)

    return curdata.drop(dropinds, axis=0)


seqinfo = pandas.read_csv(infofile, delimiter='\t', header=2)

imseq = input('Which image sequence to analyze?')

try:
    seqind = list(seqinfo[seqinfo.Sequence == imseq].index)
    if len(seqind) > 1:
        print('Multiple sequences match given name:')
        print(seqinfo[['Sequence', 'Part', 'MyFile']])
        FirstOrAll = input('Use first sequence [F], or merge all [A]?')
    else:
        FirstOrAll = 'F'

except IndexError:
    print('Sequence name does not match user input')
    print(seqinfo.Sequence)
else:
    if (FirstOrAll == 'A') | (FirstOrAll == 'F'):
        curdata = seqprocess(
                    infodf=seqinfo, ind=seqind[0], folder=parentdir,
                    scale=MicronsPerPixel)
        # Identify method for tracking blobs in moving frames.
        trackmethod = seqinfo.TrackMethod[seqind[0]]

        if FirstOrAll == 'A':
            for k in seqind[1:]:
                curdata = curdata.append(
                                seqprocess(
                                    infodf=seqinfo, ind=k, folder=parentdir,
                                    scale=MicronsPerPixel), ignore_index=True)
                if trackmethod != seqinfo.TrackMethod[k]:
                    raise SystemExit(
                        '"TrackMethod" differs among parts of image sequence.')
    else:
        raise SystemExit('Invalid choice.')

"""
CURRENTLY SET TO IGNORE DIFFERENCE BETWEEN AUTO AND MANUAL.
"""

if (trackmethod == 'Auto') | (trackmethod == 'Manual'):
    curdata['ImGroup'] = curdata['Media'].astype(int)
elif trackmethod == 'Cannot':
    # meastimelist : used to create dict of imtimedict and then associated dict
    # values with rows in curdata data frame.
    meastimelist = curdata['Time'].tolist()
    # imtimedict : used to create imgrpdict
    imtimedict = dict(zip(meastimelist,
                          curdata[['Media', 'Moving']].values))
    # imtimelist : used to create imgrpdict
    imtimelist = sorted(imtimedict.keys())
    # imgrpdict : associates image times with group numbers; group numbers
    # uniquely identify segments of the sequence in which embryos should be
    # trackable (not moving much, and in same media)
    imgrpdict = dict.fromkeys(imtimelist, 0)
    for k in range(1, len(imtimelist)):
        imgrpdict[imtimelist[k]] = imgrpdict[imtimelist[k-1]] + int(
            (imtimedict[imtimelist[k]][1] | imtimedict[imtimelist[k-1]][1]) | (
                imtimedict[imtimelist[k]][0] ^ imtimedict[imtimelist[k-1]][0]))
    # create column of curdata associating each image time with a group.
    curdata['ImGroup'] = [imgrpdict[item] for item in meastimelist]
    # cleanup
    del meastimelist, imtimelist, imgrpdict, imtimedict
else:
    raise SystemExit('Invalid track method option')

TrackPoints.linkpoints(curdata, DataColumns=['Time', 'X', 'Y', 'Major'],
               GroupNameColumn='ImGroup', BlobNameColumn='blobID', name1=0)

# # Print XY and volume data for grouped blobs.
# curgroups = curdata.groupby('blobID')
# for name, group in curgroups:
#     print(name)
#     print(group[['Media', 'Time', 'X', 'Y', 'Volume']])

## plot all linked blobs by x-y center.
#fig, ax = plt.subplots()
#colorlist = 'rgbcmyk'
#markerlist = 'o^sx+D'
#for k in set(curdata['blobID'].values.tolist()):
#    dfsub = curdata[curdata['blobID']==k]
#    colval = colorlist[k%len(colorlist)]
#    markval = markerlist[k%len(markerlist)]
#    ax.scatter(dfsub['X'].values, dfsub['Y'].values, color=colval,
#               marker=markval, alpha=0.4)

# plot all linked blobs by time-volume.
fig, ax = plt.subplots()
colorlist = 'rgbcmyk'
markerlist = 'o^sx+D'
for k in set(curdata['blobID'].values.tolist()):
    dfsub = curdata[curdata['blobID']==k]
    colval = colorlist[k%len(colorlist)]
    markval = markerlist[k%len(markerlist)]
    ax.scatter(dfsub['Time'].values, dfsub['Volume'].values, color=colval,
               marker=markval, alpha=0.4)
