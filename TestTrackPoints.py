# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 11:47:00 2016

@author: Michelangelo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy

from TrackPoints import pointcollection


#time1 = curdata[curdata['ImGroup'] == 1].Time.values[0]
#time2 = curdata[curdata['ImGroup'] == 2].Time.values[0]
#df1 = curdata[curdata['Time'] == time1]
#df2 = curdata[curdata['Time'] == time2]
#print(df1[['X', 'Y']])
#print(df2[['X', 'Y']])
#mypoints = pointcollection(
#            df1, datacols=['Time', 'X', 'Y', 'Major'], firstpointname=0)
#newpoints = pointcollection(
#            df2, datacols=['Time', 'X', 'Y', 'Major'], firstpointname=0)
#print(mypoints.matchpoints(newpoints))
#
#yspoints = copy.deepcopy(mypoints)
#yspoints.update(newpoints)

n = 4
df1 = pd.DataFrame(
                np.array([[0]*n, [0]*n, list(range(0,n)), [1]*n]).transpose(),
                columns=['T', 'X', 'Y', 'D'])
df2 = copy.deepcopy(df1)
df2['T'] = df1['T'].values + 1
df2.index = df1.index + n

mypoints = pointcollection(
            df1, datacols=['T', 'X', 'Y', 'D'], firstpointname=0)
newpoints = pointcollection(
            df2, datacols=['T', 'X', 'Y', 'D'], firstpointname=0)
yspoints = copy.deepcopy(mypoints)
yspoints.update(newpoints)

list(yspoints.points['originalInds']) == df2.index.tolist()
list(yspoints.points['names']) == list(range(0, n))
list(yspoints.points['frame']) == list(df2['T'].values)
list(yspoints.points['xs']) == list(newpoints.points['xs'])
list(yspoints.points['ys']) == list(newpoints.points['ys'])
list(yspoints.points['diams']) == list(newpoints.points['diams'])

