# -*- coding: utf-8 -*-
"""
Miscelanious functions

Created on Sun Dec 18 12:12:00 2016

@author: Michelangelo
"""
import json
import os


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
            print(savefilename +
                  ' already exists in this folder: NOT SAVED!')
    except OSError:
        with open(savefilename, 'w') as myfile:
            json.dump(mydict, myfile)
        try:
            with open(savefilename, 'r') as myfile:
                print(savefilename + ' saved.')
        except OSError:
            print('File saving error.')


def mynewfolder(mydirname, mynewfolder, maxk=10):
    """
    Create a new directory in path 'mydirname. with name 'mynewfolder'. If
    folder with name 'mynewfolder' already exists in 'mydirname', it asks to
    create a new one, incrementing from 0 to maxk; stops if it reaches maxk.

    Parameters:
    ----------
    mydirname : string
        valid path name
    mynewfolder : string
        valid folder name
    maxk : int

    Returns:
    ---------
    new folder name after incrementing (empty list - [] - if cannot create
    folder.)
    """
    if os.path.exists(os.path.join(mydirname, mynewfolder)):
        print(mynewfolder + ' already exists in this directory.')
        for k in range(maxk):
            response = input(
                    'Make new folder ' + mynewfolder + str(k) + '? y/n')
            if response is 'y' or response is 'Y':
                if os.path.exists(os.path.join(mydirname, mynewfolder +
                                               str(k))):
                    print(mynewfolder + str(k) +
                          ' already exists in this directory.')
                else:
                    mynewfolder = mynewfolder + str(k)
                    os.mkdir(os.path.join(mydirname, mynewfolder))
                    break

            elif response is 'n' or response is 'N':
                mynewfolder = []
                break
            else:
                print('Invalid response, treated as n')
                mynewfolder = []
                break

        else:
            mynewfolder = []
            print('Failed to create folder with unique name in ' +
                  str(maxk) + ' tries')

    else:
        os.mkdir(os.path.join(mydirname, mynewfolder))

    return mynewfolder
