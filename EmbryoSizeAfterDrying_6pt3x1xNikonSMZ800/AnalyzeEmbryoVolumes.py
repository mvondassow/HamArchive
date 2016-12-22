# -*- coding: utf-8 -*-
"""
Steps for analyzing embryo volume files (THIS SCRIPT DOES STEPS 4-5)
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

Created on Sun Dec 18 15:10:04 2016

@author: Michelangelo
"""
import pandas
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import time

import os, sys
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import hambits.stats as hs

from patsy import dmatrices

datafile = 'ZygoteVolumes_AveragedByRibbon.csv'

# Import selected data
consdata = pandas.read_csv(datafile)

tstyle = '%Y %B %d - %I:%M %p'
print('Run on: ', time.strftime(tstyle))
print('data file updated: ', 
      time.strftime(tstyle, time.gmtime(os.path.getmtime(datafile))))

print('data: ')
print(consdata[['File', 'TmntCat', 'Salinity', 'Volume', 'SE_Volume'
                ]].round(decimals={'Salinity': 0, 'Volume': 0, 'SE_Volume': 0,
                                 'InvSalinity': 3}))

# Generate plots
tmntgrps = consdata.groupby('TmntCat')
# Plot groups (by 1/salinity)
fig1, ax1 = plt.subplots()
for name, group in tmntgrps:
    ax1.plot(group.InvSalinity, group.Volume, marker='o', linestyle='', ms=12,
             label=name)
ax1.legend(loc='lower right')
ax1.set_ylabel('Volume, µm^3')
ax1.set_xlabel('1/salinity, (1/ppt)')
ax1.set_ylim((0, 350000))
ax1.set_xlim(0, 0.04)

# Plot groups (by salinity)
fig2, ax2 = plt.subplots()
for name, group in tmntgrps:
    ax2.plot(group.Salinity, group.Volume, marker='o', linestyle='', ms=12,
             label=name)
ax2.legend(loc='lower right')
ax2.set_ylabel('Volume, µm^3')
ax2.set_xlabel('salinity, ppt')
ax2.set_ylim((0, 350000))
ax2.set_xlim(0, 200);

# Analyze with linear regression
y, X = dmatrices('Volume ~ TmntCat + InvSalinity + InvSalinity:TmntCat',
                 data=consdata, return_type='dataframe')
mod = sm.OLS(y, X)  # Describe model
regression_result = mod.fit()  # Fit model
print(regression_result.summary())

# Check for outliers (first term in dict is number of outliers)
outlierresults = hs.GeneralizedESD(regression_result.resid.values, 10,
                                   Alpha=0.05)
print('Number of outliers detected: ', outlierresults[0])
print(pandas.DataFrame(outlierresults[1]))

# Test if makes a difference if drop value with exceptionally high salinity
consdata2 = consdata.drop(consdata[consdata['File'] == 'Z14_meas.xls'].index)
y2, X2 = dmatrices('Volume ~ TmntCat + InvSalinity + InvSalinity:TmntCat',
                   data=consdata2, return_type='dataframe')
mod2 = sm.OLS(y2, X2)  # Describe model
regression_result2 = mod2.fit()  # Fit model
print(regression_result2.summary())

# Plot with lines that were fit to all data.
xlims = [0, 0.04]
ax1.plot(xlims, regression_result.predict(np.array([[1, 0, xlims[0], 0],
                                                    [1, 0, xlims[1], 0]]
                                                    )), color='b')
ax1.plot(xlims, regression_result.predict(np.array([[1, 1, xlims[0], xlims[0]],
                                                    [1, 1, xlims[1], xlims[1]]]
                                                    )), color='g')

xlims = [0, 0.04]
# # remove last two lines added
# del ax1.lines[-2:]
ax1.plot(xlims, regression_result2.predict(np.array([[1, 0, xlims[0], 0],
                                                     [1, 0, xlims[1], 0]])),
                                                     color='b')
ax1.plot(xlims, regression_result2.predict(np.array([[1, 1, xlims[0], xlims[0]],
                                                     [1, 1, xlims[1], xlims[1]]]
                                                     )), color='g')
