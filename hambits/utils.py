# -*- coding: utf-8 -*-
"""
Miscelanious functions

Created on Sun Dec 18 12:12:00 2016

@author: Michelangelo
"""
import json


def savemydf(mydf, savefilename, extension):
    """
    Function to save Pandas data fram 'mydf' to tab-separated CSV file.
    This uses the pandas.dataframe method for saving to csv file, but wraps
    it in checks for whether the file exists and whether it's successfully
    saved.

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
