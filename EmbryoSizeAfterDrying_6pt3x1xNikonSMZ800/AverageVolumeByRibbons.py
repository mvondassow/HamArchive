# -*- coding: utf-8 -*-
"""
Steps for analyzing embryo volume files (THIS SCRIPT DOES STEPS 1-3)
0) MANUALLY IN EXCEL: Create csv file with embryo image name, treatment
(isolated embryos in salt solution or embryos in egg masses), and [expected]
salinity.
1) Read individual files:
Check that 1) have 10 lines per and correct columns: raise error if they don't.
2) Calculate XY centers of each selection and compare to selection (since index
may go off, use positional index rather than ImageJ index). Check that XY
centers are within <1 radius of each other.
3) Calculate volume and mean volume per file.
4) Plot volume vs 1/salinity
5) Linear least-squares fit to 1/salinity for each treatment (linear regression
model in statsmodel package).
6) Do fits for individual treatments and compare slopes & following Zar (1999)
to check statsmodel outputs.

Created on Thu Dec 15 20:31:54 2016

@author: Michelangelo
"""
import pandas
import numpy as np

import os, sys
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import hambits.utils as hu

umperpix = 0.548307783  # See BrightFieldVsObliqueImages.xlsx

mydir = 'Data'
mysubdir = 'Zygotes_SizeAfterDrying'
filefromlog = 'ZygoteConsolidatedInfo.csv'

consdatafile = 'ZygoteVolumes_AveragedByRibbon'

consdata = pandas.read_csv(filefromlog)


consdata = consdata[['Embryo number', 'TmntCat', 'Salinity', 'File']]
consdata.set_index('File', inplace=True, verify_integrity=True)
consdata['Volume'] = np.nan
consdata['SE_Volume'] = np.nan

# Number of rows to expect per file
numrows = 10

for curfile in consdata.index.values:
    curdata = pandas.read_csv(os.path.join(mydir, mysubdir, curfile), sep='\t')
    # Check that expected number of entries per
    if len(curdata) != numrows:
        print(curfile)
        raise SystemExit('Current file has unexpected number of rows')

    # Calculate embryo radii (measurements are in pairs: max diameter and min
    # diameter of each embryo, so Length from every other row, and divide by 4)
    radii_pix = (curdata['Length'].iloc[0::2].values +
                 curdata['Length'].iloc[1::2].values)/4
    # Check that entries are not skipped:  so the x-y centers of each ROI
    # should be within 1 embryo radius.
    curdata['Xcent'] = curdata['BX']+curdata['Width']/2
    curdata['Ycent'] = curdata['BY']+curdata['Height']/2
    dx=curdata['Xcent'].iloc[0::2].values - curdata['Xcent'].iloc[1::2].values
    dy=curdata['Xcent'].iloc[0::2].values - curdata['Xcent'].iloc[1::2].values
    checkcent = (dx*dx+dy*dy)**0.5 < radii_pix
    if not all(checkcent):
        print(curfile, checkcent)
        raise SystemExit('Measurements not in pairs?')

    # Calculate volume
    radii_um = umperpix*radii_pix
    volumes_um3 = (4*3.14159/3)*(radii_um**3)
    # Enter mean volume into consdata dataframe
    consdata.loc[curfile, 'Volume']=np.mean(volumes_um3)
    # Enter standard errors of mean volumes
    consdata.loc[curfile, 'SE_Volume']=np.std(volumes_um3
                                             )/(volumes_um3.size**0.5)

consdata['InvSalinity'] = 1/consdata['Salinity']

hu.savemydf(consdata, consdatafile, 'csv')


