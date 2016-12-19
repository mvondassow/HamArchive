# -*- coding: utf-8 -*-
"""
Various statistical tests.

Created on Sun Dec 18 11:25:56 2016

@author: Michelangelo
"""
import numpy as np
import scipy.stats as st
import warnings as wrn


def GeneralizedESD(MyData, MaxNumOutliers, Alpha=0.05):
    """
    Generalized ESD test for outliers, following:
        NIST/SEMATECH e-Handbook of Statistical Methods, "Generalized ESD Test
        for Outliers"
        http://itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm, 2016Dec17

    Parameters :
    ------------
    MyData : array-like
        MyData will be flattened if not 1D
    MaxNumOutliers : int
        upper-bound on number of suspected outliers
    Alpha : float
        Alpha level for test: desired probability of type 1 error, i.e. of
        rejecting the null hypothesis that data points are not outliers, if
        the null hypothesis (and assumptions) are true; 0 < Alpha < 1

    Returns :
    -----------
    dict :
        inds : array
            indices of outliers in MyData
        values : dict
            MaxRsInd : array
                index of maximum of standardized residuals of data in MyData,
                after removing data points corresponding to previous iterations
            vals : array
                value of array at current MaxRsInd
            MaxRs : array
                maximum of standardized residuals
            L : array
                critical value for that iteration (at Alpha)

    Example :
    ------------
    rosnerdata = [-0.25, 0.68, 0.94, 1.15, 1.20, 1.26, 1.26, 1.34, 1.38, 1.43,
                  1.49, 1.49, 1.55, 1.56, 1.58, 1.65, 1.69, 1.70, 1.76, 1.77,
                  1.81, 1.91, 1.94, 1.96, 1.99, 2.06, 2.09, 2.10, 2.14, 2.15,
                  2.23, 2.24, 2.26, 2.35, 2.37, 2.40, 2.47, 2.54, 2.62, 2.64,
                  2.90, 2.92, 2.92, 2.93, 3.21, 3.26, 3.30, 3.59, 3.68, 4.30,
                  4.64, 5.34, 5.42, 6.01]
    result = GeneralizedESD(MyData=rosnerdata, MaxNumOutliers=10, Alpha=0.05)
    print('Indices of suspected outliers: ', result[0])
    print('Table of values from calculations:')
    print(pandas.DataFrame.from_dict(result[1]))

    The output from this is (or should be) close, but not exactly the same as
    in the example from the Handbook. The difference is likely to be due to
    rounding error.

    Quote results of this example from Handbook:
        '''
              H0:  there are no outliers in the data
              Ha:  there are up to 10 outliers in the data

              Significance level:  α = 0.05
              Critical region:  Reject H0 if Ri > critical value

              Summary Table for Two-Tailed Test
              ---------------------------------------
                    Exact           Test     Critical
                Number of      Statistic    Value, λi
              Outliers, i      Value, Ri          5 %
              ---------------------------------------
                      1          3.118          3.158
                      2          2.942          3.151
                      3          3.179          3.143 *
                      4          2.810          3.136
                      5          2.815          3.128
                      6          2.848          3.120
                      7          2.279          3.111
                      8          2.310          3.103
                      9          2.101          3.094
                     10          2.067          3.085
        For the generalized ESD test above, there are essentially 10
        separate tests being performed. For this example, the largest
        number of outliers for which the test statistic is greater than
        the critical value (at the 5 % level) is three. We therefore
        conclude that there are three outliers in this data set.
        '''
    """
    MyData = np.array(MyData).astype(float)
    if MyData.ndim > 1:
        wrn.warn('Array is has ndim > 0. Results for flattened array.')
        MyData = MyData.flatten()
    n = len(MyData)

    # Initialize dict to store for values of calculations at each iteration
    IterVals={'MaxRsInd': np.full(MaxNumOutliers, np.nan),
              'MaxRs': np.full(MaxNumOutliers, np.nan),
              'Vals': np.full(MaxNumOutliers, np.nan),
              'Rcrit': np.full(MaxNumOutliers, np.nan)}
    for k in range(0, MaxNumOutliers):
        # Calculate standardized residuals (Rs), and get maximum and its index
        Rs = abs(MyData-np.nanmean(MyData))/np.nanstd(MyData)
        MaxRsInd = np.nanargmax(Rs)
        MaxRs = Rs[MaxRsInd]

        # Calculate critical value
        j = k + 1
        # Adjust 1-Alpha for 2 tails & n-j+1 possible comparisons
        Padj = 1 - Alpha/(2*(n-j+1))
        # t-value w/ cummulative probability Padj & n-1-j degrees of freedom
        Tcrit = st.t.ppf(Padj, n-j-1)
        # Calculate crtical values for standardized residual
        Rcrit = (n - j) * Tcrit / (((n - j - 1 + Tcrit**2)*(n-j+1))**0.5)

        # Store values in IterVals
        IterVals['MaxRsInd'][k] = MaxRsInd
        IterVals['MaxRs'][k] = MaxRs
        IterVals['Vals'][k] = MyData[MaxRsInd]
        IterVals['Rcrit'][k] = Rcrit

        # remove value at index for next round
        MyData[MaxRsInd] = np.nan

    # Identify indices for significant values
    OutlierInds = None
    for k in range(MaxNumOutliers-1, -1, -1):
        if IterVals['MaxRs'][k]>IterVals['Rcrit'][k]:
            OutlierInds = IterVals['MaxRsInd'][0:k+1]
            IterVals['sig'] = np.array([True]*(k+1) +
                                       [False]*(MaxNumOutliers-k-1))
            break

    return OutlierInds, IterVals


def cit(myarray, interval=0.95):
    """
    confidence interval for mean on myarray given t-distribution

##    Parameters :
    ------------


##    Returns :
    ---------

##    Example
    ----------

    """
    # Check that equivalent to one dim array
    if myarray.ndim == 1:
        # remove nans
        myarray2 = myarray[np.logical_not(np.isnan(myarray))]
        # Calculate average
        avg = np.mean(myarray2)
        # Calculate confidence interval
        lb, ub = st.t.interval(interval, len(myarray2)-1,
                               loc=avg, scale=st.sem(myarray2))
        return {'mean': avg, 'LB': lb, 'UB': ub}
    else:
        raise SystemExit('myarray should be 1 dimensional')


def cib(myarray, interval=0.95, quantile=0.5):
    """
    confidence interval for quantile (default is median) on myarray given
    binomial distribution, based on Conover, Practical Nonparametric Statistics
    1999, chapter 3.2 pp143-144

##    Parameters :
    ------------


##    Returns :
    ---------

    Example
    ----------
    From Conover p144-145
    conover = np.array([46.9, 56.8, 63.3,67.1,47.2,59.2,63.4,67.7,49.1,
                        59.9, 63.7, 73.3,56.5,63.2,64.1,78.5])
    print(cib(conover, 0.95, quantile=0.75))

##    Expected output: _____________
    """
    # Check that equivalent to one dim array
    if myarray.ndim == 1:
        # remove nans and sort
        myarray2 = np.sort(myarray[np.logical_not(np.isnan(myarray))])
        # Calculate confidence interval
        lbind, ubind = st.binom.interval(interval, len(myarray2)-1, quantile)
        lb = myarray2[int(lbind)]
        ub = myarray2[int(ubind)]
        return {'median': np.median(myarray2), 'LB': lb, 'UB': ub}
    else:
        raise SystemExit('myarray should be 1 dimensional')
