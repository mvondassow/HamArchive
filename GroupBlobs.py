# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 19:56:12 2016

@author: Michelangelo

Group and analyze blobs (embryos) from ImageJ macro output, given information
in info file.
"""

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

# QCam saves from 0 to n-1. 'end1' gives last frame (in QCams numbering) of the
# first part of the data series (before transition to high salinity); begin2
# gives first frame of the 2nd part (after transition to high salinity).
end1 = 2
begin2 = 3 

MicronsPerPixel = 0.338774005

# import data from file 
curdata = pandas.read_csv(myfolder+myfile, delimiter='\t')

# Split out frame (image) and time information from Label column 
curdata['Image'] = pandas.to_numeric(
    curdata['Label'].str.rsplit(':',1).str[1].str.rsplit('_',1).str[0])
curdata['Time'] = pandas.to_numeric(curdata['Label'].str.rsplit('_',1).str[1])

# add columns assigning membership to first or second part of data series.
curdata['part1'] = curdata['Image'] <= end1
curdata['part2'] = curdata['Image'] >= begin2

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