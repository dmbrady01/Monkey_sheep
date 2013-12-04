function [] = db_preparepsth(window,spike_count,name)
%db_preparesmoothpsth Prepares figure for smooth PSTH
%   A series of commands to standardize the smooth PSTH and label the axis.

smooth_max = 1.05*max([max(spike_count.(name).ipsi) max(spike_count.(name).contra)]); %Calculates the max spike/sec for the ipsi and contra
%eyes to set the y axis for the subplots. Y max will be the max firing rate
%multiplied by 1.05 (can change if you want it to be smaller or larger)

%Sets the x axis to include prestimulus, stimulus, and poststimulus time
xlim([-window.prestim (window.stimulus+window.prestim)])
xlabel('Time (msec)') %label for x axis

%Sets the y axis from 0 spikes/sec to smooth_max
ylim([0 smooth_max])
ylabel('Spikes/sec') %label for y axis

%Draws lines for the beginning and end of the stimulus
%Start of stimulus line
line([0 0], [0 smooth_max],... %coordinates (time is 0)
    'Color', [0 0 0],... %color of line, currently black
    'LineWidth', 2) %width of line

line([window.stimulus window.stimulus], [0 smooth_max],... %coordinates (time is at the end of the stimulus)
    'Color', [0 0 0],... %color of line, currently black
    'LineWidth', 2,... %width of line
    'LineStyle', '--') %make it a dotted line

end

