#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
db_spikewidth.py: Python script that processes csv files of spike waveform data
to determine if a cell is fast-spiking or not.

db_spikewidth Calculates whether an average spike waveform should be considered
fast (lower than fast_threshold) or regular (above regular_threshold) based
on its peak to trough time.

One inputs a path to the folder with spike waveforms csvs. The output will be a
csv file with the cell name and its classification
"""


__author__ = "DM Brady"
__datewritten__ = "22 May 2018"
__lastmodified__ = "22 May 2018"


# Required modules
import pandas as pd
import os
import sys


# Load data from either excel or csv
def load_excel_file(file):
    """Loads excel or csv file"""
    if '.csv' in file:
        return pd.read_csv(file)
    elif '.xlsx' in file:
        return pd.read_excel(file)
    else:
        raise TypeError('You need to pass a csv or xcel file')

def check_for_waveform_data(file, term='Spk'):
    if (term in file) and (('.csv' in file) or ('.xlsx' in file)):
        return True
    else:
        return False

def walkthrough_directories(folder):
    tree = [(root, subdir, file) for root, subdir, file in os.walk(folder)]
    folders_with_data = []
    for folder in tree:
        if any(check_for_waveform_data(file) for file in folder[2]):
            folders_with_data.append(folder[0])
    return folders_with_data

# Classify waveform
def classify_waveform(dataframe, fast_threshold=400, regular_threshold=400):
    # Get average waveform
    average_waveform = dataframe.mean(axis=0)
    # Find trough
    trough = average_waveform.idxmin()
    # Find peak (after trough)
    peak = average_waveform.loc[average_waveform.index > trough].idxmax()
    # Find peak to trough duration
    peak_trough_duration = peak - trough
    # Classify
    if peak_trough_duration <= fast_threshold:
        return 'Fast'
    elif peak_trough_duration > regular_threshold:
        return 'Regular'
    else:
        return 'Cannot classify'

# Gets list of files that are spikewaveform data
def classify_folder_of_waveforms(folder, fast_threshold=400, regular_threshold=400):
    """Gets a list of waveform files."""
    data_folders = walkthrough_directories(folder)
    for data_folder in data_folders:
        contents = os.listdir(data_folder)
        results = []
        for file in contents:
            if check_for_waveform_data(file):
                df = load_excel_file(data_folder + os.sep + file)
                neuron_type = classify_waveform(df, fast_threshold, regular_threshold)
                name = file.split('_')[0]
                name = file.split('.')[0]
                results.append((name, neuron_type))
        if len(results) > 0:
            df = pd.DataFrame(results, columns=['Name', 'Type'])
            df.to_csv(data_folder + os.sep + 'spikewidth.csv')


if __name__ == '__main__':
    try:
        folder = sys.argv[2]
    except:
        folder = '/Users/DB/Development/Monkey_sheep/data/Spike_waveform_files/'
    classify_folder_of_waveforms(folder, fast_threshold=400, regular_threshold=400)

