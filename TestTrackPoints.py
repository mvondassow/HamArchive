# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 11:47:00 2016

@author: Michelangelo
"""

import pandas
import numpy as np
import copy

import TrackPoints
import importlib
import warnings

#importlib.reload(TrackPoints)
warnings.filterwarnings('error')


def test_identicaldfs(nmax=5, verbose=False):
    """
    Test if can match points between identical dataframes.

    Returns true if TrackPoints.pointcollection.update successfully makes
    pointcollection dfs equal after updating, when point positions and numbers
    are identical between input dfs; also tests that initial pointcollection df
    and final pointcollection df are equal except in the infocols.
    """
    nleft = nmax
    xy = np.random.random((nmax, 2))

    df1 = pandas.DataFrame({'T': np.zeros(nmax), 'X': xy[:, 0], 'Y': xy[:, 1],
                            'D': np.ones(nmax)}, index=np.arange(0, nmax))
    df2 = pandas.DataFrame({'T': np.ones(nleft), 'X': xy[:nleft, 0],
                            'Y': xy[:nleft, 1], 'D': np.ones(nleft)},
                           index=np.arange(nmax, nmax+nleft))

    mypoints = TrackPoints.pointcollection(df1, datacols=['X', 'Y'],
                                           infocols=['T', 'D'],
                                           firstpointname=0, weights=[1, 1])
    newpoints = TrackPoints.pointcollection(df2, datacols=['X', 'Y'],
                                            infocols=['T', 'D'],
                                            firstpointname=0, weights=[1, 1])

    yspoints = copy.deepcopy(mypoints)
    yspoints.update(newpoints)

    out1 = np.all(newpoints.points.values == yspoints.points.values)

    if verbose:
        print(df1)
        print(df2)
        print(mypoints.points)
        print(yspoints.points)

    return out1


def test_lostsome(nstart=5, nmissing=1, verbose=False):
    """
    Test if can find points when points are missing, but xy positions do not
    change between dataframes df1 and df2.

    Returns true if TrackPoints.pointcollection.update successfully makes
    pointcollection dfs equal after updating, except for having one extra point
    in the updated df (yspoints.points); also tests that initial
    pointcollection df and final pointcollection df are equal except in the
    infocols.
    """
    nleft = nstart-nmissing
    xy = np.random.random((nstart, 2))

    df1 = pandas.DataFrame({'T': np.zeros(nstart), 'X': xy[:, 0],
                            'Y': xy[:, 1], 'D': np.ones(nstart)},
                           index=np.arange(0, nstart))
    df2 = pandas.DataFrame({'T': np.ones(nleft), 'X': xy[:nleft, 0],
                            'Y': xy[:nleft, 1], 'D': np.ones(nleft)},
                           index=np.arange(nstart, nstart + nleft))

    mypoints = TrackPoints.pointcollection(df1, datacols=['X', 'Y'],
                                           infocols=['T', 'D'],
                                           firstpointname=0, weights=[1, 1])
    newpoints = TrackPoints.pointcollection(df2, datacols=['X', 'Y'],
                                            infocols=['T', 'D'],
                                            firstpointname=0, weights=[1, 1])
    yspoints = copy.deepcopy(mypoints)
    yspoints.update(newpoints)

    if verbose:
        print(df1)
        print(df2)
        print(mypoints.points)
        print(yspoints.points)

    yspoints.points.sort_values('names', axis=0, inplace=True)
    newpoints.points.sort_values('names', axis=0, inplace=True)
    mypoints.points.sort_values('names', axis=0, inplace=True)

    # Test that yspoints contains the points in newpoints
    out1 = np.all(newpoints.points.drop('names', axis=1).values ==
                  yspoints.points.drop('names', axis=1).iloc[:-nmissing,
                                                             :].values)
    # Test that all points in mypoints are preserved in yspoints,
    # including points missing from new points 
    out2 = np.all(mypoints.points.drop(
                  ['T', 'inputInds'], axis=1).values ==
                  yspoints.points.drop(['T', 'inputInds'], axis=1).values)

    # Test that names are unique
    out3 = len(yspoints.points['names'].values) == len(set(
                yspoints.points['names'].values))

    return out1, out2, out3


def test_gainedsome(nstart=5, ngained=1, verbose=False):
    """
    Test if can find points when some points areadded, but xy positions do not
    change between dataframes df1 and df2.

    Returns true if TrackPoints.pointcollection.update successfully makes
    pointcollection dfs equal after updating, except for having one extra point
    in the updated df (yspoints.points); also tests that initial
    pointcollection df and final pointcollection df are equal except in the
    infocols.
    """
    nfinal = nstart + ngained
    xy = np.random.random((nfinal, 2))

    df1 = pandas.DataFrame({'T': np.zeros(nstart), 'X': xy[:nstart, 0],
                            'Y': xy[:nstart, 1], 'D': np.ones(nstart)},
                            index=np.arange(0, nstart))
    df2 = pandas.DataFrame({'T': np.ones(nfinal), 'X': xy[:, 0],
                            'Y': xy[:, 1], 'D': np.ones(nfinal)},
                            index=np.arange(nstart, nstart + nfinal))

    mypoints = TrackPoints.pointcollection(df1, datacols=['X', 'Y'],
                                           infocols=['T', 'D'],
                                           firstpointname=0, weights=[1, 1])
    newpoints = TrackPoints.pointcollection(df2, datacols=['X', 'Y'],
                                            infocols=['T', 'D'],
                                            firstpointname=0, weights=[1, 1])
    yspoints = copy.deepcopy(mypoints)
    yspoints.update(newpoints)

    if verbose:
        print(df1)
        print(df2)
        print(mypoints.points)
        print(yspoints.points)

    yspoints.points.set_index('X', inplace=True)
    newpoints.points.set_index('X', inplace=True)
    mypoints.points.set_index('X', inplace=True)

    # Test that yspoint contains all points from newpoints.
    out1 = all([np.all(newpoints.points.drop('names', axis=1).loc[
                x, :].values == yspoints.points.drop('names', axis=1
                ).loc[x, :].values) for x in newpoints.points.index.values])
    # Test that yspoints contains all points from mypoints
    out2 = all([np.all(mypoints.points.drop(['T', 'inputInds'], axis=1
                       ).loc[x, :].values == yspoints.points.drop(
                       ['T', 'inputInds'], axis=1).loc[x, :].values)
                       for x in mypoints.points.index.values])
    # Test that names of points are unique
    out3 = len(yspoints.points['names'].values) == len(set(
                yspoints.points['names'].values))

    return out1, out2, out3


def test_reordered(nstart=5, ngained=1, verbose=False):
    """
    Test if can find points when points are reordered and some may be added,
    but xy positions do not change between dataframes df1 and df2.

    Returns true if TrackPoints.pointcollection.update successfully makes
    pointcollection dfs equal after updating, except for having extra points
    in the updated df (yspoints.points); also tests that initial
    pointcollection df and final pointcollection df are equal except in the
    infocols.
    """
    nfinal = nstart + ngained
    xy = np.random.random((nfinal, 2))

    df1 = pandas.DataFrame({'T': np.zeros(nstart), 'X': xy[:nstart, 0],
                            'Y': xy[:nstart, 1], 'D': np.ones(nstart)},
                            index=np.arange(0, nstart))
    df2 = pandas.DataFrame({'T': np.ones(nfinal), 'X': xy[:, 0],
                            'Y': xy[:, 1], 'D': np.ones(nfinal)},
                            index=np.arange(nstart, nstart + nfinal))
    df2.sort_values('Y', inplace=True)
    df2.reset_index(drop=True, inplace=True)

    mypoints = TrackPoints.pointcollection(df1, datacols=['X', 'Y'],
                                           infocols=['T', 'D'],
                                           firstpointname=0, weights=[1, 1])
    newpoints = TrackPoints.pointcollection(df2, datacols=['X', 'Y'],
                                            infocols=['T', 'D'],
                                            firstpointname=0, weights=[1, 1])
    yspoints = copy.deepcopy(mypoints)
    yspoints.update(newpoints)

    if verbose:
        print(df1)
        print(df2)
        print(mypoints.points)
        print(yspoints.points)

    yspoints.points.set_index('X', inplace=True)
    newpoints.points.set_index('X', inplace=True)
    mypoints.points.set_index('X', inplace=True)

    # Test that yspoint contains all points from newpoints.
    out1 = all([np.all(newpoints.points.drop('names', axis=1).loc[
                x, :].values == yspoints.points.drop('names', axis=1
                ).loc[x, :].values) for x in newpoints.points.index.values])
    # Test that yspoints contains all points from mypoints
    out2 = all([np.all(mypoints.points.drop(['T', 'inputInds'], axis=1
                       ).loc[x, :].values == yspoints.points.drop(
                       ['T', 'inputInds'], axis=1).loc[x, :].values)
                       for x in mypoints.points.index.values])
    # Test that names of points are unique
    out3 = len(yspoints.points['names'].values) == len(set(
                yspoints.points['names'].values))

    return out1, out2, out3


def test_varypositions(nstart=3, ngained=2, relativenoise=0.1, verbose=False):
    """
    NOTE: THIS TEST WILL NOT NECESSARILY YIELD TRUE WHEN THE FUNCTION WORKS
    PROPERLY BECAUSE THE FUNCTION SEEKS TO MIMIMIZE DISTANCES BETWEEN OLD & NEW
    POINTS. THE JIGGLE IN POSITION MAY ALTER WHICH POINTS ARE CLOSEST.

    Test if can find points when point positions are jiggled, indices are
    reordered, and some points are added.

    Returns true if TrackPoints.pointcollection.update successfully makes
    pointcollection dfs equal after updating, except for having one extra point
    in the updated df (yspoints.points); also tests that initial
    pointcollection df and final pointcollection df are equal except in the
    infocols.
    """
    nfinal = nstart + ngained
    xy = np.random.random((nfinal, 2))
    xy2 = xy + relativenoise * np.random.random((nfinal, 2))

    df1 = pandas.DataFrame({'T': np.zeros(nstart), 'X': xy[:nstart, 0],
                            'Y': xy[:nstart, 1], 'D': np.arange(0, nstart)},
                            index=np.arange(0, nstart))
    df2 = pandas.DataFrame({'T': np.ones(nfinal), 'X': xy2[:, 0],
                            'Y': xy2[:, 1], 'D': np.arange(0, nfinal)},
                            index=np.arange(nstart, nstart + nfinal))
    df2.sort_values('Y', inplace=True)
    df2.reset_index(drop=True, inplace=True)

    # This time use column D as a persistant index to use when testing matching
    mypoints = TrackPoints.pointcollection(df1, datacols=['X', 'Y'],
                                           infocols=['T', 'D'],
                                           firstpointname=0, weights=[1, 1])
    newpoints = TrackPoints.pointcollection(df2, datacols=['X', 'Y'],
                                            infocols=['T', 'D'],
                                            firstpointname=0, weights=[1, 1])
    yspoints = copy.deepcopy(mypoints)
    yspoints.update(newpoints)

    if verbose:
        print(df1)
        print(df2)
        print(mypoints.points)
        print(yspoints.points)

    
    yspoints.points.set_index('D', inplace=True)
    newpoints.points.set_index('D', inplace=True)
    mypoints.points.set_index('D', inplace=True)

    # Test that yspoint names contains all points from newpoints.
    out1 = all([np.all(newpoints.points.drop('names', axis=1).loc[
                d, :].values == yspoints.points.drop('names', axis=1
                ).loc[d, :].values) for d in newpoints.points.index.values])
    # Test that yspoints contains all points from mypoints: ONLY MATCH 'names'
    # COLUMN
    out2 = all([mypoints.points['names'].loc[d] == 
                yspoints.points['names'].loc[d]
                for d in mypoints.points.index.values])
    # Test that names of points are unique
    out3 = len(yspoints.points['names'].values) == len(set(
                yspoints.points['names'].values))

    return out1, out2, out3


test_varypositions(nstart=3, ngained=2, relativenoise=0.1, verbose=True)
