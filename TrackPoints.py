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
from copy import deepcopy


class pointcollection:
    def __init__(self, df, datacols=['X', 'Y'], infocols=['T', 'D'],
                 firstpointname=0, weights=[1, 1]):
        """
        Initialize point collection

        Parameters
        ----------
        df : pandas.core.frame.DataFrame
            Must have columns corresponding to datacols
        datacols : str
            Names of columns containing data to use for matching points
        infocols : str
            Names of columns containing other information to keep with points
        firstpointname : int
            Numeric name to assign first point in df
        weights : list of ints, same number of elements as datacols
            weighting to apply to columns for matching points when doing sum
            of squares (e.g. 1*(deltaX**2)+1*(deltaY**2))

        pointcollection.points : pandas dataframe with columns datacols and
            infocols (see above), 'names' (numeric names for points), and
            'inputInds' (index from dataframe used
            to create pointcollection; this column gets updated by
            pointcollection.update(newpointcollection) to match inputInds in
            newpointcollection for all points that match between the two
            pointcollections).
        pointcollection.cols : dict of strings
            Contains names of datacolumns and infocolumns
        pointcollection.weights = list of weights to apply to datacolumns when
            calculating sum of squares
        pointcollection.nextpointname : int
            starting point for naming new points added to
            pointcollection.points upon pointcollection.update.
        """
        self.cols = {"datacolumns": datacols, "infocolumns": infocols}
        self.weights = weights
        self.nextpointname = firstpointname + len(df.values)

        self.points = deepcopy(df[datacols + infocols])
        self.points['inputInds'] = np.array(df.index)
        self.points['names'] = np.arange(firstpointname, self.nextpointname)

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
            rowinds/colinds : row/col indices in matrix of shape datashape
            that correspond to indices in list q of flattened array.
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
        for p in r:
            indlist += [p*datashape[1] + k for k in range(
                                                        0, datashape[1])]
        return indlist

    @classmethod
    def matchfun(cls, initialmatch, rankinds_in, datashape):
        """
        Pass through matching next nearest pairs of points (i.e. lowest ranks
        of squared distances).

        Parameters
        ----------
        cls : class
        initialmatch : int
            index for element in rankinds_in corresponding to first matched
            pair of points
        rankinds_in : list of ints
            contains ranks of elements in sqrdists (array of squared distances)
            for all pairs of points from self.points and newpoints.points.
        datashape : tuple
            shape (#rows, #columns) of original array of squared distances
            (sqrdists: number of rows = number of points in self.points;
            number of columns = number of points in newpoints.points)

        Returns :
        dict with keys:
        'matched' : indices in rankinds_in corresponding to matched points
        'lost' : indices of points that appear in self.points but not in
            newpoints.points
        'new' : indices of points that appear in newpoints.points but not in
            self.points
        'conflicts' : indices of point pairs in which one point cannot be
            matched to its nearest neighbor because its nearest neighbor is
            matched to a different point (e.g. a is closer to b than any other
            point, but b is closer to c than any point).
        """
        # Use this to copy values, not references.
        rankinds_out = [q for q in rankinds_in]
        match = [initialmatch]
        conflicts = []  # Initialize array to store conflicting point linkings
        iterations = 0  # Use to prevent looping too long.
        keepgoing = True
        while (len(rankinds_out) > 0) & (iterations <= datashape[1]
                                         ) & keepgoing:
            iterations += 1
            # Get rid of elements in rankinds corresponding to same column as
            # match[0]
#            print('A', rankinds_out, match, datashape, iterations)
            matchcolinds = cls.FlatIndsInSameColumn(match[-1], datashape)
            for q in matchcolinds:
#                print('B', matchcolinds, q, datashape, rankinds_out, match)
                try:
                    rankinds_out.remove(q)
                except:
                    print("Couldn't find q to remove. How did I get here?")
                    pass

            # Find next match
            conflicts = []  # store conflicts where col ind is missing
            for k in rankinds_out:
#                print('k=', k)
                if k in cls.FlatIndsInSameRow(match, datashape):
                    conflicts += [k]
                    # Create stop condition: goes false if conflict for last
                    # item in rankinds_out
                    keepgoing = k != rankinds_out[-1]
#                    print('got to conflicts: ', conflicts)
                else:
                    match += [k]
#                    print('got to match')
                    break

        return {"matched": match, "conflicts": conflicts}

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
            keys: matched, new, lost, conflicts, sqrdists
            For all keys except sqrdists, values are tuples of two lists. The
                first list contains indices for new points; the 2nd list
                contains indices for old points (for new and lost, one entry
                will be []).
            'matched' : indices in corresponding to matched points
            'lost' : indices of points that appear in self.points but not in
                newpoints.points
            'new' : indices of points that appear in newpoints.points but not
                in self.points
            'conflicts' : indices of point pairs in which one point cannot be
                matched to its nearest neighbor because its nearest neighbor is
                matched to a different point (e.g. a is closer to b than any
                other point, but b is closer to c than any point).
            'sqrdists' : squared distances between centers of each pair of
                matched points.
        """
        # Calculate array of squared distances btwn pointcollection &
        # newpoints. For large arrays, it might be better not to calculate all
        # elements, but for arrays of expected size it should be faster.
        if self.points.columns.tolist() != newpoints.points.columns.tolist():
            raise ValueError("Old & new pointcollection columns don't match")

        olddata = (self.points[self.cols['datacolumns']].values)*np.array(
                                                        [self.weights])
        newdata = (newpoints.points[self.cols['datacolumns']].values)*np.array(
                                                        [self.weights])
        for c in range(0, olddata.shape[1]):
            # Using numpy's None indexing to expand array dimensions to make
            # a column vector from olddata and a row vector from newdata. Numpy
            # propagates to form rxc array upon addition.
            if c == 0:
                sqrdists = (olddata[:, c, None]-newdata[None, :, c])**2
            else:
                sqrdists += (olddata[:, c, None]-newdata[None, :, c])**2
#        print(sqrdists)

        rankinds = sqrdists.argsort(axis=None).tolist()

        # Starting point of matching: link pair of points with shortest
        # distance between pointcollection and newpoints
        match = rankinds[0]

        # Find first iteration of matches among points.
        matchdict = self.matchfun(match, rankinds, sqrdists.shape)

        # If any conflicts (more than one point has same nearest neighbor),
        # test if alternative choices of starting point give lower sum of
        # squared distances ('cumsum').
        if len(matchdict['conflicts']) > 0:
            cursum = sqrdists.flat[matchdict["matched"]].sum()
            for conf in matchdict['conflicts']:
                newmatchdict = self.matchfun(
                    conf, rankinds, sqrdists.shape)
                newsum = sqrdists.flat[newmatchdict["matched"]].sum()
                if newsum < cursum:
                    cursum = newsum
                    matchdict['matched'] = newmatchdict['matched']

        # Convert indices in squared distance matrix to row-col indices: row
        # indices (R) are indices of old points; col indices (C) refer to new
        # points.
        matchRCs = self.FlatIndsToRowCol(matchdict['matched'], sqrdists.shape)
        newCs = ([], list(set(
                range(0, sqrdists.shape[1])).difference(set(matchRCs[1]))))
        lostRs = (list(set(
                range(0, sqrdists.shape[0])).difference(set(matchRCs[0]))), [])
        conflictRCs = self.FlatIndsToRowCol(matchdict['conflicts'],
                                            sqrdists.shape)

        return {'matched': matchRCs, 'new': newCs, 'lost': lostRs,
                'conflicts': conflictRCs, 'sqrdists': sqrdists[matchRCs]}

    def update(self, newpoints, verbose=False):
        """
        Modify values in self.points to those from newpoints,points where
        points match, and add points to self.points that appear in newpoints
        but not in self.points.

        newpoints : pointcollection object
        verbose : boolean
        """
        # Find point matches.
        matchdict = self.matchpoints(newpoints)
        if verbose:
            print(matchdict)

        # deepcopy newpoints.points to avoid changing original
        newptsdfcopy = deepcopy(newpoints.points)

        # Get position index of 'names' column; raise error if columns do not
        # match between self and new pointcollections
        ptcols = self.points.columns.tolist()
        if ptcols != newptsdfcopy.columns.tolist():
            raise ValueError("Old & new pointcollection columns don't match")

        namepos = ptcols.index('names')
        # For matched points, change name of points in newptsdfcopy. Can't do
        # all at once without difficult to identify error (probably due to
        # modifying df as it goes through the list.
        for spts, npts in zip(
                            matchdict['matched'][0], matchdict['matched'][1]):
            newptsdfcopy.iloc[npts, namepos] = self.points.iloc[spts, namepos]

        # Rename new points in newptsdfcopy and update nextpointname
        for npt in matchdict['new'][1]:
            newptsdfcopy.iloc[npt, namepos] = self.nextpointname
            self.nextpointname += 1

        # Use concatenate to add rows from lost points in old points df.
        self.points = pandas.concat((
                    self.points.iloc[matchdict['lost'][0], :], newptsdfcopy))


def linkpoints(df, DataColumns=['X', 'Y'], InfoColumns=['Time', 'Major'],
               GroupNameColumn='ImGroup', BlobNameColumn='blobID', name1=0,
               ColWeights=[1, 1]):
    """
    Link points in data frame df
    """
    if type(name1) != int:
        YorN = input(
             'name1 must be an int. Continue with firstpointname=0? y/n')
        if YorN == 'n':
            raise ValueError('First point name not an int')
        else:
            name1 = 0

    newevent = True

    gnc = GroupNameColumn  # Column containing group names for images
    fc = InfoColumns[0]  # Column containing frame info (e.g. time of shot)

    bnc = BlobNameColumn  # Column to contain names for points/blobs
    # Initialize column for point/blob names: using None made assignment crash
    df[bnc] = -float('inf')

    p1 = name1  # Initialize firs point/blob name
    for imgroup in sorted(set(df[gnc].values)):
        f1 = df[df[gnc] == imgroup].loc[:, fc].values[0]
        dfinit = df[df[fc] == f1]
        initialpoints = pointcollection(dfinit, datacols=DataColumns,
                                        infocols=InfoColumns,
                                        firstpointname=p1, weights=ColWeights)
        # Turns out Python sets aren't sorted even though they print in order,
        # so have to sort here to go through images in order.
        for fnew in sorted(set(df[df[gnc] == imgroup].loc[:, fc])):
            dfnew = df[df[fc] == fnew]
            newpoints = pointcollection(dfnew, datacols=DataColumns,
                                        infocols=InfoColumns,
                                        firstpointname=p1, weights=ColWeights)
            initialpoints.update(newpoints)
            for k in df[df[fc] == fnew].index:
                #try:
                df.loc[k, bnc] = initialpoints.points[initialpoints.points[
                                        'inputInds'] == k].names.values
#                except ValueError:
#                    if newevent:
#                        newevent = False
#                        print(fc, '==', fnew, '; index=', k)
#                        print("Points =", initialpoints.points)
#                        print('New points=', newpoints.points)
#                        print("df IDs =", df.loc[k, bnc])
#                        break

        p1 = max(initialpoints.points["names"]) + 1

    return df
