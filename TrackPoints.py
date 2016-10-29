# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 17:01:33 2016

Class and function to link points by position between two dataframes.

I've written it for the specific application tracking Ham. embryos by x-y
coordinates while keeping track of frame time and embryo diameter, but
might be helpful to rewrite for more general columns, and potentially used
more than the x-y column similarity to match points (e.g. weighted sum of
squares of another quantifier (e.g. diameter)).

Storing values in pointcollection.points as dataframe might also make indexing
easier...

@author: Michelangelo
"""

import pandas
import numpy as np


class pointcollection:
    def __init__(self, dataframe, datacols=['Frame', 'X', 'Y',
                                            'Diameter'], firstpointname=0):
        """
        Initialize point collection

        Parameters
        ----------
        dataframe : pandas.core.frame.DataFrame
            Must have columns corresponding to datacols
        datacols : str
            Names of columns containing frame names or times, coordinates and
            diameters of points
        firstpointname : int
            Numeric name to assign first point in dataframe

        Returns
        -------
        numpy.ndarray : contains names (ints) for points
        """
        cols = {'frame': datacols[0], 'x': datacols[1], 'y': datacols[2],
                'diams': datacols[3]}

        # Values in dict self.points should iterables of the same length
        self.points = {'names': np.arange(firstpointname,
                                       firstpointname + len(dataframe.values)),
                    'frame': dataframe[cols['frame']].values,
                    'xs': np.expand_dims(dataframe[cols['x']].values, axis=1),
                    'ys': np.expand_dims(dataframe[cols['y']].values, axis=1),
                    'diams': dataframe[cols['diams']].values,
                    'originalInds': np.array(dataframe.index)}
        self.cols = cols

    @staticmethod
    def FlatIndsToRowCol(q, datashape):
        """
        Converts indices in flattened 2D array into (row, col) indices of
        2D array.

        Parameters :
        ------------
        q : int or iterable of ints
        datashape : tuple

        returns :
        ---------
        tuple (length 2)
        """
        if not hasattr(q, '__iter__'):
            q = [q]
        rowinds = [item // datashape[1] for item in q]
        colinds = [item % datashape[1] for item in q]
        return rowinds, colinds

    @staticmethod
    def RowColToFlatInds(rowcol, datashape):
        """
        Convert (row, col) index tuple to index for flattend array.
        2D ARRAYS ONLY!

        Parameters
        ----------
        rowcol : list, tuple, or array; numeric
            row number in position 0, column in position 1
        datashape : tuple, with two elements.
            element 0: number of rows, element 1: number of columns in
            array
        """
        return rowcol[0]*datashape[1]+rowcol[1]

    @classmethod
    def FlatIndsInSameColumn(cls, q, datashape):
        """
        Given 'q', an index value from a flat array, return the other
        indices in a 2D array of shape 'datashape' that would be in the
        same column as q.
        """
        r, c = cls.FlatIndsToRowCol(q, datashape)
        indlist = []
        for p in c:
            indlist += [p + datashape[1] * k for k in range(
                                                        0, datashape[0])]
        return indlist

    @classmethod
    def FlatIndsInSameRow(cls, q, datashape):
        """
        Given 'q', an index value from a flat array, return the other
        indices in a 2D array of shape 'datashape' that would be in the
        same column as q.
        """
        r, c = cls.FlatIndsToRowCol(q, datashape)
        indlist = []
        for q in r:
            indlist += [q*datashape[1] + k for k in range(
                                                        0, datashape[1])]
        return indlist

    @classmethod
    def matchfun(cls, initialmatch, rankinds_in, datashape):
        """
        Pass through matching next nearest pairs of points (i.e. lowest ranks of
        squared distances).
        """
        # Use this to copy values, not references.
        rankinds_out = [q for q in rankinds_in]
        match = [initialmatch]
        conflicts = []  # Initialize array to store conflicting point linkings
        iterations = 0  # Use to prevent looping too long.
        while (len(rankinds_out) > 0) & (
                                    iterations <= (datashape[0]*datashape[1])):
            iterations += 1
            # Get rid of elements in rankinds corresponding to same column as
            # match[0]
#           print('A', rankinds, match, datashape)
            matchcolinds = cls.FlatIndsInSameColumn(match[-1], datashape)
            for q in matchcolinds:
#                print('B', matchcolinds, q, datashape, rankinds, match)
                try:
                    rankinds_out.remove(q)
                except:
                    pass

            # Find next match
            conflicts = []  # store conflicts where col ind is missing
            for k in rankinds_out:
                if k in cls.FlatIndsInSameRow(match[0], datashape):
                    conflicts += [k]
                else:
                    match += [k]
                    break

            # Break when run out of possibilities for matching
            if len(match) < iterations:
                break

        return match, rankinds_out, conflicts

    def matchpoints(self, newpoints):
        """
        Match points described in dataframe to points in point collection.
        Dataframe must matches format of dataframe used to generate
        pointcollection.

        Uses greedy algorithm to link points in pointcollection and newpoints
        The way I am implementing this ignores possiblities of ties (they
        will be rare in this application) but they could be dealt with either
        by checking for them, or by repeating with noise added, or ...

        Parameters
        -----------
        newpoints : pointcollection

        Returns
        -------
        dict :
        keys: matched, unmatched, conflicts, sqrdists
        matched/unmatched: indices of matched/unmatched points; conflicts:
        indices where two new points share same nearest neighbor in old points.
        Each is a tuple of two lists. The first list contains indices for old
        points; the 2nd list contains indices for new points. 'sqrdists':
        squared distances between centers of each pair of matched points.
        """
        # Calculate array of squared distances btwn pointcollection &
        # newpoints. For large arrays, it might be better not to calculate all
        # elements, but for arrays of expected size it should be faster.
        sqrdists = (self.points['xs'] - newpoints.points['xs'].transpose()
            )**2 + (self.points['ys'] - newpoints.points['ys'].transpose())**2
#        print(sqrdists)

        rankinds = sqrdists.argsort(axis=None).tolist()

        # Starting point of matching: link pair of points with shortest distance
        # between pointcollection and newpoints
        match = rankinds[0]

        # Find first iteration of matches among points.
        match, unmatched, conflicts = self.matchfun(
            match, rankinds, sqrdists.shape)

        # If any conflicts (more than one point has same nearest neighbor), test
        # if alternative choice of starting point gets lower sum of squared
        # distances ('cumsum').
#        print('Match ', match, ' Unmatch ', unmatched, ' Conflicts ', conflicts)
        if len(conflicts) > 0:
            cursum = sqrdists.flat[match].sum()
#            print('Match sum ', cursum)
            for k in conflicts:
                newmatch, newunmatched = self.matchfun(
                    conflicts[k], rankinds, sqrdists.shape)[0:2]
#                print(self.matchfun(conflicts[k], rankinds, sqrdists.shape))
#                print(rankinds)
                newsum = sqrdists.flat[newmatch].sum()
#                print('Newmatch sum ', newsum, ' Newmatch ', newmatch, ' New unmatched ', newunmatched)
                if newsum < cursum:
                    cursum = newsum
                    match, unmatched = newmatch, newunmatched

#        print('Final sum ', cursum, ' Final match ', match, ' Final unmatched ', unmatched)
        # Convert indices in squared distance matrix to row-col indices: row
        # indices are indices of old points; col indices refer to new points.
        matchrcs = self.FlatIndsToRowCol(match, sqrdists.shape)
        unmatchrcs = self.FlatIndsToRowCol(unmatched, sqrdists.shape)

        return {'matched': matchrcs, 'unmatched': unmatchrcs,
                'conflicts': conflicts, 'sqrdists': sqrdists.flat[match]}

    def update(self, newpoints):
        matchdict = self.matchpoints(newpoints)

        for key in self.points.keys():
            if key != 'names':
                self.points[key][matchdict['matched'][0]] = newpoints.points[
                    key][matchdict['matched'][1]]
                self.points[key] = np.concatenate((self.points[key],
                            newpoints.points[key][matchdict['unmatched'][1]]))
            else:
                self.points['names'] = np.concatenate((self.points['names'],
                                    np.arange(self.points['names'][-1]+1,
                                          self.points['names'][-1] + 1 + len(
                                          matchdict['unmatched'][1])) ))


def linkpoints(df, DataColumns=['Time', 'X', 'Y', 'Major'],
               GroupNameColumn='ImGroup', BlobNameColumn='blobID', name1=0):
    """
    Link points in data frame df
    """
    if type(name1) != int:
        YorN = input(
             'name1 must be an aint. Continue with firstpointname=0? y/n')
        if YorN == 'n':
            raise ValueError('First point name not an int')
        else:
            name1 = 0

    newevent = True

    gnc = GroupNameColumn  # Column containing group names for images
    fc = DataColumns[0]  # Column containing frame info (e.g. time of shot)

    bnc = BlobNameColumn  # Column to contain names for points/blobs
    df[bnc] = None  # Initialize column for point/blob names

    p1 = name1  # Initialize firs point/blob name
    for imgroup in set(df[gnc].values):
        f1 = df[df[gnc] == imgroup].loc[:, fc].values[0]
        dfinit = df[df[fc] == f1]
        initialpoints = pointcollection(dfinit, datacols=DataColumns,
                                        firstpointname=p1)
        for fnew in set(df[df[gnc] == imgroup].loc[:, fc]):
            dfnew = df[df[fc] == fnew]
            newpoints = pointcollection(dfnew, datacols=DataColumns,
                                        firstpointname=p1)
            initialpoints.update(newpoints)
            for k in df[df[fc] == fnew].index:
                try:
                    df.loc[k, bnc] = initialpoints.points["names"][
                            initialpoints.points[
                            "originalInds"].tolist().index(k)]
                except ValueError:
                    if newevent:
                        newevent = False
                        print(fc, '=', fnew, '; index=', k)
                        print("Points =", initialpoints.points)
                        print('New points=', newpoints.points)
                        print("df IDs =", df.loc[k, bnc])
                        break

        p1 = max(initialpoints.points["names"]) + 1

    return df

