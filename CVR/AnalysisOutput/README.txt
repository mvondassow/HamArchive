CVR output files:

-----------

CVR_cleavers_summaryinfo.csv and CVR_Zygotes_summaryinfo.csv
Comma separated values files from data frame with columns :
    FileName
    TimeConstUpper : upper bound on time constant
    VolRatio : min volume/max volume
    RecoveredFraction : 
    MinDelay : time between last image from first image group (in original media) to time of first image in final image group: represents minimum measurable delay time)
    TimeToCutoff : Time after media change at which volume crosses cutoff (cutoff volume = min volume + (initial volume - min volume)*1-e^(-ntcs); ntcs = number of time constants; fOR ntcs = 2 -- as used for files created on 2016Dec27 -- this corresponds to the volume falling >/=86.5% of the way to the minimum volume).

Due to a quirk in the Pandas package, the columns are not ordered in the obvious way.

-----------

CleaverSummary_CIs.json and ZygoteSummary_CIs.json:
json formatted files (openable as text) with summary info, including parameters used in the analysis, and CIs (UB: lower bound; UB: upper bound) on the mean or median (method of CI calculation shown).

Parameters: 'tou': conversion of time values in column 'Time' to desired units (minutes/second), 'ipu': images per time unit (e.g frames per minute), 'tmax': how many units forward to calculate for recovered volume fraction (60 min), 'ntcs': number of time constants to use when calculating bound on time constant for shrinkage (2 time constants corresponds to 1-exp(-2) ~ 0.865: for ntcs=2 the script 'CalculationsForCVR.py' identifies when the volume falls by 86.5% or more of the way to the minimum volume).

Due to a quirk in how Python handles dictionaries, the order of items in curly braces varies unpredictably.