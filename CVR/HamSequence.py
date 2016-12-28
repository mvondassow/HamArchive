# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 19:56:12 2016

@author: Michelangelo

'''
2016 Dec 26: Discovered and corrected two errors in volume calculation: forgot
to divide average diameter by 2 (i.e. sum_of_diameters/4), and missed a zero
in conversion to nL for figure. Now converted to pL instead.

Need to update to re-run.
'''

Group and analyze blobs (embryos) from ImageJ macro output, given information
in info file. Put into dataframe 'curdata'; process, link blobs, and plot.

Classes and functions:
---------------------------
BlobViewer class & methods :
    Takes dataframe (curdata) and folder information, reads in image files, and
    creates a figure window that user can advance/decrement through to see
    whether blobs/ROIs were linked correctly.
seqprocess
    Import and process data from ImageJ macro.

Requires:
--------
•Tab delimited file (with appropriate format and labels) containing info needed
to identify breaks in sequence, or measurements to treat specially.
•Tab delimited file with appropriate columns and labels containing data for
each sequence
•Scale info (CURRENTLY CODED IN SCRIPT)
•Directory path for files (CURRENTLY CODED IN SCRIPT)
Modules/packages:
    pandas, numpy, json, matplotlib.pyplot, skimage.external.tifffile,
    TrackPoints

Steps to process files from ImageJ results
------------------------------------------
1) Read data file and info file.

2) Calculations on data and split or remove data based on info file.
 2a: Split 'Label' to get time stamps and file codes
 2b: Calculate volume
 2c: Divide data into two sets (before and after transition) and delete frames
 during transition (identify manually)
 2d: For sequences split into separate parts in different polders, combine data
 from separate parts indo dataframe (curdata) if requested.
 2d: Divide images into groups based on tracking method user entered in info
 file

3) Match ROIs (blobs) from one frame to next based on overlap in XY center. In
 at least one case, an embryo was missed in one frame (crossed edge) so join if
 overlap is close to blob that appeared earlier in same image group.

4) Plot volume data.

5) Show image sequence with labeled blobs
"""
import pandas
import numpy as np
import json
import TrackPoints
import matplotlib.pyplot as plt
import skimage.external.tifffile as tifmod

import os, sys
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import hambits.utils as hu

# Parent directory
parentdir = 'C:\\Users\\Michelangelo\\Documents\\Ham\\'
# Name of file with info about sequences
infofile = os.path.join(parentdir, 'CellVolumeRegulation.txt')

# Scale in images
MicronsPerPixel = 0.338774005


class BlobViewer:
    def __init__(self, blobdata, folderlist, scale=10):
        """
        Create a event-driven figure with which one can scroll through images
        in image sequence and see blobs labeled by 'blobID' at XY centers.
        Image files are found in folders of folderlist, and referred to by info
        in blobdata dataframe.

        Actions (when mouse is within figure axes) :
        -------
        '.' key press : show next image in sequence, & corresponding blob IDs
        ',' key press : show previous image in sequence, & its blob IDs
        'c' key press : close figure window

        Parameters :
        ------------
        blobdata : pandas dataframe
            Must contain columns 'X', 'Y', 'Label' (with the form '*:*:*'),
            'Time', and 'blobID'; blobID, X, & Y must be numeric.
        folderlist :  list
            List of folder paths (as strings) associated with image sequences
            described in blobdata
        scale : int
            How much to downsample the image by.

        Returns :
        ---------
        None
        """
        self.foldernames = folderlist
        self.blobdata = blobdata
        self.scale = int(scale)

        # Timelist contains times associated with images, in order.
        self.timelist = sorted(list(set(blobdata['Time'].values)))
        self.t = self.timelist[0]

        # Create figure window and initialize axis (just to simplify code)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        # Connect user actions (key presses) with a function (self.wherenext)
        self.useraction = self.fig.canvas.mpl_connect(
            'key_press_event', self.wherenext)

        # Initialize dicts to store file names (self.filedict), indices in
        # blobdata df for blobs in each image (self.indsdict), and image data
        # (self.imdict) as nparrays.
        self.filedict = dict.fromkeys(self.timelist)
        self.indsdict = dict.fromkeys(self.timelist)
        self.imdict = dict.fromkeys(self.timelist)

        # Go through times in timelist and get indices for blobs, filename
        # and image
        for loopt in self.timelist:
            self.indsdict[loopt] = self.blobdata[self.blobdata['Time'] == loopt
                                                 ].index.tolist()
            self.filedict[loopt] = self.blobdata.loc[min(self.indsdict[loopt]),
                                                     'Label'].rsplit(
                                                     ':', 1)[1] + '.tif'
            # Try to find filename for image at time loopt in each folder; stop
            # when find it.
            for myfolder in self.foldernames:
                try:
                    currentimage = tifmod.imread(os.path.join(myfolder,
                                                 self.filedict[loopt]))
                except:
                    continue
                else:
                    self.imdict[loopt] = currentimage[::self.scale,
                                                      ::self.scale]
                    break

        # Create figure with first image and associated blobs.
        self.showblobs()

    def showblobs(self):
        """
        Update figure with image and blobs associated with time self.t

        Parameters :
        ------------
        self : BlobViewer object

        Returns :
        ---------
        None
        """
        self.ax.remove()
        self.ax = self.fig.add_subplot(111)
        self.ax.imshow(self.imdict[self.t], cmap='gray')
        self.ax.set_title('Image: ' + self.filedict[self.t])
        for ind in self.indsdict[self.t]:
            self.ax.text(self.blobdata.loc[ind, 'X']/self.scale,
                         self.blobdata.loc[ind, 'Y']/self.scale,
                         int(self.blobdata.loc[ind, 'blobID']),
                         color=(1, 1, 0))
        plt.draw()

    def wherenext(self, event):
        """
        Update self.t and show blobs and image associated with new self.t
        based on user action (using connection created in __init__).
        """
        tind = self.timelist.index(self.t)
        if event.inaxes:
            if event.key == '.':
                try:
                    self.t = self.timelist[tind + 1]
                except:
                    self.t = self.timelist[0]
                self.showblobs()
            elif event.key == ',':
                self.t = self.timelist[tind - 1]
                self.showblobs()
            elif event.key == 'c':
                plt.close(self.fig)
            else:
                pass


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
    seqfolder = os.path.join(folder, infodf.SetDir[ind],
                             infodf.SequenceDir[ind], 'flattened')
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
    curdata = pandas.read_csv(os.path.join(seqfolder, seqfile), delimiter='\t')

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
        scale*(curdata['Feret']+curdata['MinFeret'])/4)**3

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

def savecurdata(datadf, infodf, ind):
    if len(ind) > 0:
        ind = ind[0]
    newfile = infodf.MyFile[ind]

"""
Read datafile with user-generated metadata about each image sequence and get
user input for which sequence to process
"""
seqinfo = pandas.read_csv(infofile, delimiter='\t', header=2)
imseq = input('Which image sequence to analyze?')

"""
Read in files associated with the seqence 'imseq'. For any sequence which are
split into different parts, ask if should combine the parts.
The get the trackmethod (whether to try to track blobs in images listed in
infofile as 'moving'.
Then process (and possibly combine) information from csv file containing info
on blobs.
"""
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
Group data into image groups (in colum 'ImGroup' in curdata df). 
If 'trackmethod' set to 'Auto' or 'Manual', just groups by media; if
'trackmethod' set to 'cannot', this creates separate groups for different media
and for each images marked as moving, so that no blob will be linked between
groups.
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

"""
Link blobs in image groups (column 'ImGroup') using TrackPoints module
"""
TrackPoints.linkpoints(curdata, DataColumns=['X', 'Y'],
                       InfoColumns=['Time', 'Major'],
                       GroupNameColumn='ImGroup', BlobNameColumn='blobID',
                       name1=0, ColWeights=[1, 1])

"""
Plot all linked blobs by time-volume.
"""
fig, ax = plt.subplots()
colorlist = 'rgbcmyk'
markerlist = 'o^sx+D'
for k in set(curdata['blobID'].values.tolist()):
    dfsub = curdata[curdata['blobID'] == k]
    colval = colorlist[int(k) % len(colorlist)]
    markval = markerlist[int(k) % len(markerlist)]
    # Time of media change was between 0 and 1 minute after last frame in first
    # media, therefore assign time of media change to midpoint (+30s), although
    # the actual media change may have been a bit faster (not accounting for
    # mixing time).
    StartTime = curdata[curdata['Media'] == 0].Time.iloc[-1] + 30
    ax.scatter((dfsub['Time'].values-StartTime)/60,
               dfsub['Volume'].values/(1000),
               color=colval, marker=markval, alpha=0.4)
ax.set_ylim(bottom=0, top=600)
ax.set_xlim(left=-10, right=70)
ax.set_ylabel('Volume, pL')
ax.set_xlabel('Time after media change, min')
ax.set_title(imseq + ': cell volume over time')


"""
Create figure to check link among blobs.
"""
foldernames = [parentdir + seqinfo.loc[q, 'SetDir'] + '\\' +
               seqinfo.loc[q, 'SequenceDir'] + '\\flattened\\'
               for q in seqind]
temp = BlobViewer(curdata, foldernames)

"""
Function to save curdata to tab separated CSV file.
Too much of a fight to get matplotlib to display figures before going on to the
rest of the script and fucking up.
When ready, run: hu.savemydf(curdata, destfilename, destextension)
"""
destfilename = seqinfo.loc[seqind, 'MyFile'].values[0].split('.')[0] + \
               '_Processed'
destextension = 'csv'
