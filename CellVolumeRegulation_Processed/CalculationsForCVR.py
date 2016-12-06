# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 18:33:54 2016
Calculate parameters from cell volume regulation experiment:

For each, use mean values per image because embryos are not the same before and
after media change.
1) Biological question: Does envelope or cell membrane provide permeability
barrier at time scale of exposure?
Empirical question: What is ime scale of volume loss? Is it smaller that the
time scale of exposure in tides (>=1 hr)?
1a) Calculate ratio of volume lost at 3 min to maximum volume lost (CI based on
normal distribution or sign test/mann-whitney U?). Is majority of volume lost
by 3 min? [SKIP: harder to interpret than 1b because arbitrary time scale]
1b) Calculate upper bound on time constant: time point when volume lost >= 0.63
of max volume lost: get CI (assume log-normal distribution because can't be <0;
or based on sign-test/mann-whitney U: least restrictive assumptions). Does CI
overlap 1 hr?
    2016Dec02: Before calculating stats, I realized it would be better to
    calculate the time constant upper bound as 1/2 of the time at which:
    volume_lost>=(1-exp(-2))*max_volume_lost (e.g. volume lost corresponds to
    two time constants) this would provide the narrowest bounds on the time
    constant.
    2016Dec02: CI based on binomial distribution following Conover, Practical
    Nonparametric Statistics 1999

2) Biological question: Is it as expected if the embryo behaves as a simple
osmometer with 100% osmotically active fraction?
Empirical question: is the minimum relative volume comparable to the ratio of
salinities? (Perfect osmometer with 100% osmotically active, expect:
    Smedia2/Smedia1~Vembryo2/Vembryo1; Smedia1/Smedia2 = 1.5)
Relative volume: min volume/initial volume; CI based on t-test or
mann-whitney-U? Does CI overlap Smedia1/Smedia2?

3) Biological question: Could cell volume regulation compensate for water loss
at relevant time scale of exposure (>=1hr)?
Empirical question: What is an upper bound on volume regulation?
Calculate: volume recovered/max volume lost;

volume recovered = vol lost at 1 hr - max vol lost;
vol lost = initial vol - vol at time t
initial vol = average of vol at first 3 time points.


Rib02: Noticable blip in calculated volume at ~image 15 because a tiny bit of
debris got lumped into ROI... Maybe an argument for using fit elipse when
calculating volumes instead of Feret's diameters?

@author: Michelangelo
"""
import numpy as np
import pandas
import scipy.stats as st
import time, json

ZygoteFiles = {'rib01': 'Results_CVR_rib01_MeasEmbV5_centroid_Processed.csv',
               'rib02': 'Results_CVR_rib02_MeasEmbV5_centroid_Processed.csv',
               'rib03': 'Results_CVR_rib03_MeasEmbV5_centroid_Processed.csv',
               'rib04': 'Results_CVR_rib04_MeasEmbV5_centroid_Processed.csv',
               'rib05': 'Results_CVR_rib05_MeasEmbV5_masked_centroid_Processed.csv',
               'rib06': 'Results_CVR_rib06_MeasEmbV5_centroid_Processed.csv'}
CleaverFiles = {'rib07': 'Results_CVR_rib07_MeasEmbV5_centroid_Processed.csv',
                'rib08': 'Results_CVR_rib08part1_MeasEmbV5_centroid_Processed.csv',
                'rib09': 'Results_CVR_rib09_MeasEmbV5_masked_centroid_Processed.csv',
                'rib10': 'Results_CVR_rib10_MeasEmbV5_centroid_Processed.csv',
                'rib11': 'Results_CVR_rib11_MeasEmbV5_centroid_Processed.csv',
                'rib12': 'Results_CVR_rib12_MeasEmbV5_masked_centroid_Processed.csv',
                'rib13': 'Results_CVR_rib13_MeasEmbV5_centroid_Processed.csv'}

ZygoteTransitions = {'rib01': 2, 'rib02': 2, 'rib03': 2, 'rib04': 2,
                     'rib05': 2, 'rib06': 2}
CleaverTransitions = {'rib07': 2, 'rib08': 2, 'rib09': 2, 'rib10': 2,
                      'rib11': 2, 'rib12': 2, 'rib13': 2}

# Dict of constants for time conversions: 'tou' time values in column 'Time' to
# desired units, ipu : images per time unit, 'tmax' : how many units forward to
# calculate, 'ntcs' : number of time constants to use when calculating bound on
# time constant.
myparams = {'tou': 1/60, 'ipu': 1, 'tmax': 60, 'ntcs': 2}


def cvrsummarydata(stagefiles, stagetransitions, params):
    """
    Calculate summary values for all files named in 'stage'

    Parameters :
    ------------
    stagefiles : dict
        each value is name of file; each file must be in working directory and
        be a csv with columns including 'Time', 'Image', 'ImGroup', 'blobID',
        'Media', and 'Volume'; all contain numeric datatypes.
    stagetransitions : dict
        each value is the last image name before the treatment/media transition
    timeparams : dict
        Dict of constants for time conversions: 'tou' time values in column
        'Time' to desired units, ipu : images per time unit, 'tmax' : how many
        units forward to calculate, 'ntcs' : number of time constants.

    Returns
    -------
    Data frame with columns :
        FileName
        TimeConstUpper : upper bound on time constant
        VolRatio : min volume/max volume
        RecoveredFraction : (vol(60s) - min(vol(t)))/(initial vol - min(vol(t))
    """
    filekeys = sorted(stagefiles.keys())
    nkeys = len(filekeys)
    for k in range(0, nkeys):
        summarydict = summarize(stagefiles[filekeys[k]],
                                stagetransitions[filekeys[k]],
                                params)
        if k == 0:
            summarydf = pandas.DataFrame(index=range(0, nkeys),
                                         columns=summarydict.keys())
        for key in sorted(summarydict.keys()):
            summarydf.loc[k, key] = summarydict[key]
                        

    return summarydf


def summarize(curfile, transitionimage, params):
    """
    Calculate summary values for file named curfile

    Parameters :
    ------------
    curfile : string
        name of a csv file in the working directory with columns including
        'Time', 'Image', 'ImGroup', 'Media', and 'Volume'; all with numeric
        datatypes.
    stagetransitions : dict
        each value is the last image name before the treatment/media transition
    params : dict
        Dict of constants for conversions: 'tou' time values in column
        'Time' to desired units, ipu : images per time unit, 'tmax' : how many
        units forward to calculate; 'ntcs' : number of time constants

    Returns
    -------
    Data frame with columns :
        FileName
        TimeConstUpper : upper bound on time constant
        VolRatio : min volume/max volume
        RecoveredFraction : (vol(60s) - min(vol(t)))/(initial vol - min(vol(t))
        MinDelay : time of first
    """
    metadatacols = ['Image', 'Media', 'ImGroup']

    # Assumes all columns are numeric.
    curdata = pandas.read_csv(curfile)
    curdata['Time'] = curdata['Time'].values*params['tou']
    grpd = curdata[metadatacols+['Time', 'Volume']].groupby(curdata['Time'])
    # Check that only one value per group for metadata
    if np.any(grpd[metadatacols].max().values-grpd[metadatacols].min().values):
        print(curfile)
        print(curdata[metadatacols])
        raise SystemExit('Metadata column(s) with multiple values in a group')
    else:
        grpd = grpd.mean()
        mediavals = set(grpd['Media'].values)
        if len(mediavals) != 2:
            print(curfile)
            print("Media values: ", mediavals)
            raise SystemExit("Wrong number of values in column 'Media'")
        else:
            # Split df into set in first treatment, and set in last image group
            # (want to use only last image group because comparison seems
            # better if use same embryos within treatment 2).
            initialdf = grpd[grpd['Media'] == min(mediavals)]
            finaldf = grpd[grpd['ImGroup'] == max(grpd['ImGroup'].values)]

            # Calculate time when media changed.
            ttransition = float(initialdf[initialdf['Image'] == transitionimage
                                          ].index.values)

            # Find minimum volume and time of minimum volume
            volumes = finaldf['Volume'].values
            times = finaldf['Time'].values
            indmin = np.argmin(volumes)
            minvol = volumes[indmin]
            tofmin = times[indmin] - ttransition

            # Initial volume (for embryos in first media/treatment)
            initialvol = initialdf['Volume'].mean()
            # Maximum volume lost
            vollost = initialvol - minvol

            # Calculate upper bound on time constant for volume loss.
            # Cutoff volume for estimation of upper bound on time constant.
            # transitionimage is last image in media 0. tcross is first image
            # in which embryos have lost >=0.63 of maximum volume lost.
            cutoffvol = initialvol - vollost*(1-np.exp(params['ntcs']))
            # tcross is time when volume reaches cutoff volume.
            tcross = min(times[volumes < cutoffvol])

            # tcub is upperbound on time constant
            tcub = (tcross - ttransition)/params['ntcs']

            # Calculate minimum relative volume
            minrelvol = minvol/initialvol

            # Calculate fraction of volume regained.
            try:
                tend = ttransition + params['tmax']
                tendind = findnearest(times, tend,
                                      params['ipu']/2)
                # find volume at tendind
                endvol = volumes[tendind]
                recfraction = (endvol - minvol)/vollost
            except:
                recfraction = np.nan

            # Calculate time between first usable frame in second media &
            # transition time
            mindelay = min(times)-ttransition
            return {'FileName': curfile, 'MinVolRatio': minrelvol,
                    'TimeConstEst': tcub, 'RecoveredFraction': recfraction,
                    'TimeOfMinVol': tofmin, 'MinDelay': mindelay}


def findnearest(myarray, target, bounds):
    """
    Find index of element closest to target, that is also within +/- bounds of
    target.
    This function may be slow for sorted arrays.
    """
    newarray = abs(myarray - target)
    candidate = np.argmin(newarray)
    if newarray[candidate] < bounds:
        return candidate
    else:
        return np.nan


def cit(myarray, interval=0.95):
    """
    confidence interval for mean on myarray given t-distribution
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

## Test example from Conover p144-145
#conover = np.array([46.9, 56.8, 63.3,67.1,47.2,59.2,63.4,67.7,49.1,59.9, 63.7, 
#                    73.3,56.5,63.2,64.1,78.5])
#print(cib(conover, 0.95, quantile=0.75))

def savemydf(mydf, savefilename, extension):
    """
    Function to save mydf to tab separated CSV file.
    
    Parameters
    ----------
    mydf : dataframe
    savefilename : string, name of file without extension
    extension : string, name of extension
    """
    saved = False
    while not saved:
        try:
            with open(savefilename + '.' + extension, 'r') as myfile:
                pass
            print(savefilename + ' already exists in this folder: NOT SAVED!')
            again = input('Try again using new name? y/n')
            if again == 'y':
                savefilename = input('Enter new file name (w/o extension).')
            elif again == 'n':
                saved = True
                break
            else:
                print('invalid entry. Try again.')
                continue
        except OSError:
            mydf.to_csv(savefilename + '.' + extension)
            try:
                with open(savefilename + '.' + extension, 'r') as myfile:
                    pass
                print(savefilename + '.' + extension + ' saved.')
                saved = True
                break
            except OSError:
                print('File saving error.')


def savedictasjson(mydict, savefilename):
    """
    Save ditionary as json file
    """
    try:
        with open(savefilename, 'r') as myfile:
            pass
        print(savefilename + ' already exists in this folder: NOT SAVED!')
    except OSError:
        with open(savefilename, 'w') as myfile:
            json.dump(mydict, myfile)
        try:
            with open(savefilename, 'r') as myfile:
                pass
            print(savefilename + ' saved.')
        except OSError:
            print('File saving error.')


# Generate summary data and save files
ZygoteSummary = cvrsummarydata(ZygoteFiles, ZygoteTransitions, myparams)
savemydf(ZygoteSummary, 'CVR_zygotes_summaryinfo', 'csv')
CleaverSummary = cvrsummarydata(CleaverFiles, CleaverTransitions, myparams)
savemydf(CleaverSummary, 'CVR_cleavers_summaryinfo', 'csv')

# Calculte CIs for parameters of interest and save in json format.
savefilename = 'CVR_summary_data.json'
descriptions = {'RecoveredFraction': {'about':
                                      '(V(tmax)-min(V))/(V(initial)-min(V))',
                                      'ci method': 'T', 'ci interval': 0.95},
                'MinVolRatio': {'about': 'min(V)/V(initial)',
                                'ci method': 'T', 'ci interval': 0.95},
                'TimeConstEst': {'about':
                                 'Upper bound on time constant for shrinkage',
                                 'ci method': 'Bernoulli', 'ci interval': 0.95
                                 }}

for item in (CleaverSummary, ZygoteSummary):
    if item is CleaverSummary:
        myinfo = ['CleaverSummary']
    elif item is ZygoteSummary:
        myinfo = ['ZygoteSummary']
    else:
        print('Unexpected item name')
    myfilename = myinfo[0] + '_CIs.json'
    myinfo += ['Run on: ' + time.ctime()]
    myinfo += [myparams]
    for key in descriptions:
        mydata = pandas.to_numeric(item[key].values)
        if descriptions[key]['ci method'] == 'T':
            myinfo += [{key: [descriptions[key],
                              cit(mydata, descriptions[key]['ci interval'])]}]
        elif descriptions[key]['ci method'] == 'Bernoulli':
            myinfo += [{key: [descriptions[key],
                              cib(mydata, descriptions[key]['ci interval'])]}]
        else:
            print('problem with ', key)
    savedictasjson(myinfo, myfilename)


"""
23Nov2016: Checked calculations manually for rib08. But does it make sense to
look at recovered fraction from time of transition, rather than time since min
volume?
Checked CIs in Excel
"""
