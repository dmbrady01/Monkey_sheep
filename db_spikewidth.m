function [] = db_spikewidth(path, fast_threshold=400, regular_threshold=400)
%db_spikewidth Calculates whether an average spike waveform should be considered
%fast (lower than fast_threshold) or regular (above regular_threshold) based
%on its peak to trough time.
%
%One inputs a path to the folder with spike waveforms csvs. The output will be a
%csv file with the cell name and its classification

