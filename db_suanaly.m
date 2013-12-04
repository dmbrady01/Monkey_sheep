%Single Unit Analysis
%
%Calculates baseline firing rate and evoked firing rate for ipsi and
%contralateral eyes. Calculates OD score. Makes peristimulus time
%histograms (smooth and bin).
%
%Written by DM Brady 7/12


%% What to do before using this program
%This program is designed to analyze single unit data from a multielectrode
%array in response to visual stimuli (using Visage). This program assumes
%that you recorded from the contralateral eye first. If you did not, please
%note that everything listed as contra or ipsi will be backwards.

%After sorting your cells using the online/offline sorters, open
%Neuroexplorer. Select all the single units of interest and the event1
%checkbox (that has the timing for the visual stimuli). Click on the rate
%histogram (upper right corner), make sure the data binning is 1 ms
%(0.001s), and click on the matlab tab. Make sure the outgoing data file is
%called data and click send to matlab.

%Once the file is in matlab, simply type the name of this program:
%db_suanaly



%% Make folder to save data
%This first cell is to make a folder for all of your data and figures. The
%program will ask a series of questions to make your folder and will also
%ask for the coordinates of your recording. Note that when it asks for your
%name, write it the same everytime otherwise it will make a bunch of
%folders with different versions of your name.


user = input('What is your name?  ', 's'); %asks for your name to make a folder

folder_name = input('What is the name of your experiment?  ', 's'); %asks for the name of your experiment
date_of_exp = input('When did you do this experiment? \n(example: 06Jul2012)  ', 's'); %asks for the date of your experiment

%coordinates will be a structure with three different fields corresponding
%to the A/P, M/L, and depth of your recording. Type coordinate.(field of
%interest) to look at that number.
coordinates.AP = input('What are your anterior/posterior coordinates?  '); %asks for the anterior/posterior coordinates
coordinates.ML = input('What are your medial/lateral coordinates?  '); %asks for the medial/lateral coordinates
coordinates.depth = input('What is the depth of your probe?  '); %asks for the depth of your probe

mkdir(['D:/' user '/' folder_name '_' date_of_exp]) %the actual function
%that will create the folder to store all of your information. Currently,
%it will save in the D:/ drive. Ex:
%  D:/(user)/(folder name)_(date of experiment)
%  D:/Daniel/LynxKO_05Jul2012



%% Convert single unit data from plexon into a structure of 'timestamps'
%This cell will calculate the number of visual stimuli used and convert
%data from a sequence of 0's and 1's to a structure with the times at which
%events occured (whether they are spikes or the presentation fo a visual
%stimulus).


% Trials_per_eye = 2; %This variable asks for how many times you present each stimulus to each eye. If you change the
% %visage program, make sure to change this too!
% 
% GAngle = [0;90;180;270]; %This matrix lists the orientations you used in your visual stimulus program. This information is 
% %currently not used, but you could use it to calculate orientation
% %selectivity in the future.
% 
% numberoftrials = Trials_per_eye*length(GAngle); %Calculates the total number of visual stimuli each eye received.

%The numberoftrials per eye
% numberoftrials = 8;
% timestamps = db_convert2timestamps(data,numberoftrials);

[timestamps, numberoftrials] = db_convert2timestamps2(data);%Calls the function db_convert2timestamps2 which transforms
%the variable 'data' from your recording to the times when spikes or visual
%stimuli occured.



%% Create window structure for calculating firing rates and making peristimulus histograms
%The structure window will have values for the stimulus length, pre and
%post stim period for making histograms, and the size of bins for the
%smooth and binned PSTH

%Stimulus length and pre/post stimulus period
window.stimulus = 4000; %number of milliseconds the stimulus is on
window.firing_window = window.stimulus; %window (in msecs) from the stimulus onset that we use to measure firing rate. If you want to use the entire
%time of the stimulus, set it equal to window.stimulus
window.prestim = 4000; %number of milliseconds before stimulus to calculate baseline (also used to plot post stimulus spikes)

%Smooth PSTH 
window.window = 500; %length (in msec) of bin for psth
window.step = 2; %time difference between centers of windows

%Bin PSTH
window.bin = 500; %length (in msec) of bin for histogram



%% Calculate firing rates (evoked and baseline)
%This for loop will create a structure called firingrate with many fields.
%The first field will be the neuron (labeled neuron1 to however many you
%have). For each neuron, there will be ipsi, contra, and OD score fields.
%Within the ipsi and contra fields are evoked and baseline fields. Within
%the baseline and evoked fields are the spikes/sec per trial, the average
%firing rate, and the standard deviation fo the firing rate.

eye_open = {'contra' 'ipsi'};

for i = 1:size(data,2)-1 %to repeat the following loop for the number of neurons recorded in this penetration
    
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron
    
    %loop for either ipsi or contra eye
    for k = 1:length(eye_open)
    
        for j = 1:numberoftrials
            %calculates evoked firing rate for either ipsi or contra eye for each trail by
            %finding and adding all the 1's (meaning a spike occured) from the
            %start of a visual stimulus (timestamps.events.ipsi) until the end
            %of the stimulus(window.stimulus). This is divided by the stimulus
            %length(in msecs)and then divided by 1000 to convert to spikes/sec.
            firingrate.(name).(eye_open{k}).evoked.per_trial(j) = sum(data(timestamps.events.(eye_open{k})(j):timestamps.events.(eye_open{k})(j)...
                +window.firing_window,i))./((window.firing_window)./1000);
            
            %lists when the spikes occured (msec after event)
            firingrate.(name).(eye_open{k}).evoked.per_trial_latencies{j} = find(data(timestamps.events.(eye_open{k})(j):timestamps.events.(eye_open{k})(j)...
                +window.stimulus,i));
            
            %latency to first spike
            if ~isempty(firingrate.(name).(eye_open{k}).evoked.per_trial_latencies{j})
                firingrate.(name).(eye_open{k}).evoked.latency_first_spike(j) = firingrate.(name).(eye_open{k}).evoked.per_trial_latencies{j}(1);
            else
                firingrate.(name).(eye_open{k}).evoked.latency_first_spike(j) = NaN;
            end
            
            %lists when the spikes occured including prestim time (for rasters)
            firingrate.(name).(eye_open{k}).evoked.for_rasters{j} = find(data(timestamps.events.(eye_open{k})(j)-window.prestim...
                :timestamps.events.(eye_open{k})(j)...
                +window.stimulus+window.prestim,i))-window.prestim;
            
            
            %calculates the average evoked firing rate for the ipsi or contra eye by
            %averaging the per trial firing rates from above.
            firingrate.(name).(eye_open{k}).evoked.average = mean(firingrate.(name).(eye_open{k}).evoked.per_trial);
            
            %calculates the standard deviation of the evoked firing rates for
            %the ipsi or contra eye
            firingrate.(name).(eye_open{k}).evoked.sd = std(firingrate.(name).(eye_open{k}).evoked.per_trial);
            
            %calculates baseline firing rate for ipsi or contra eye using the same process
            %as the evoked rate. The baseline period is the prestimulus period
            %specified (window.prestim).
            firingrate.(name).(eye_open{k}).baseline.per_trial(j) =...
                sum(data(timestamps.events.(eye_open{k})(j)-window.firing_window-1:timestamps.events.(eye_open{k})(j)-1,i))...
                ./(window.firing_window./1000);
%             firingrate.(name).(eye_open{k}).baseline.per_trial_latencies{j} =...
%                 find(data(timestamps.events.(eye_open{k})(j)-window.prestim-1:timestamps.events.(eye_open{k})(j)-1,i));
            firingrate.(name).(eye_open{k}).baseline.average = mean(firingrate.(name).(eye_open{k}).baseline.per_trial);
            firingrate.(name).(eye_open{k}).baseline.sd = std(firingrate.(name).(eye_open{k}).baseline.per_trial);
        end
        
        %median latency to first spike
        %get rid of NaN trials
        firingrate.(name).(eye_open{k}).evoked.latency_first_spike = firingrate.(name).(eye_open{k}).evoked.latency_first_spike(...
            ~isnan(firingrate.(name).(eye_open{k}).evoked.latency_first_spike));
        %calculates median latency
        firingrate.(name).(eye_open{k}).evoked.median_latency = median(firingrate.(name).(eye_open{k}).evoked.latency_first_spike);
        
        %Calculates whether cell fires significantly from baseline
        [firingrate.(name).(eye_open{k}).paired_test.h firingrate.(name).(eye_open{k}).paired_test.p...
            firingrate.(name).(eye_open{k}).paired_test.CI] = ttest(firingrate.(name).(eye_open{k}).evoked.per_trial,...
            firingrate.(name).(eye_open{k}).baseline.per_trial);
        
        if firingrate.(name).(eye_open{k}).paired_test.h == 1
            firingrate.(name).(eye_open{k}).sig_diff = 'yes';
        elseif firingrate.(name).(eye_open{k}).paired_test.h == 0
            firingrate.(name).(eye_open{k}).sig_diff = 'no';
        end
        
    end
    
%     %similar for loop as above but for contra eye.
%     for j = 1:length(timestamps.events.contra)       
%         %evoked firing rate for contra eye
%         firingrate.(name).contra.evoked.per_trial(j) = sum(data(timestamps.events.contra(j):timestamps.events.contra(j)+window.stimulus,i))./(window.stimulus./1000);
%         firingrate.(name).contra.evoked.per_trial_latencies{j} = find(data(timestamps.events.contra(j):timestamps.events.contra(j)+window.stimulus,i));
%         firingrate.(name).contra.evoked.average = mean(firingrate.(name).contra.evoked.per_trial);
%         firingrate.(name).contra.evoked.sd = std(firingrate.(name).contra.evoked.per_trial);
%         
%         %latency to first spike
%        if ~isempty(firingrate.(name).contra.evoked.per_trial_latencies{j})
%             firingrate.(name).contra.evoked.latency_first_spike(j) = firingrate.(name).contra.evoked.per_trial_latencies{j}(1);
%        else
%            firingrate.(name).contra.evoked.latency_first_spike(j) = NaN;
%        end
% 
%     
%         %contra baseline firing rate
%         firingrate.(name).contra.baseline.per_trial(j) = sum(data(timestamps.events.contra(j)-window.prestim-1:timestamps.events.contra(j)-1,i))./(window.prestim./1000);
%         firingrate.(name).contra.baseline.per_trial_latencies{j} = find(data(timestamps.events.contra(j)-window.prestim-1:timestamps.events.contra(j)-1,i));
%         firingrate.(name).contra.baseline.average = mean(firingrate.(name).contra.baseline.per_trial);
%         firingrate.(name).contra.baseline.sd = std(firingrate.(name).contra.baseline.per_trial);
%     end
%     
%     %median latency to first spike
%     %get rid of NaN trials
%     firingrate.(name).contra.evoked.latency_first_spike = firingrate.(name).contra.evoked.latency_first_spike(...
%         ~isnan(firingrate.(name).contra.evoked.latency_first_spike));
%     %calculates median latency
%     firingrate.(name).contra.evoked.median_latency = median(firingrate.(name).contra.evoked.latency_first_spike);
    
    %Calculate OD Score using the formula below
    % ODscore = [(ipsi evoked - ipsi baseline) - (contra evoked - contra baseline)] / [(ipsi evoked - ipsi baseline) + (contra evoked - contra baseline)]
    % Score of -1 is all contra, score of 1 is all ipsi
    firingrate.(name).ODscore = ((firingrate.(name).ipsi.evoked.average - firingrate.(name).ipsi.baseline.average)-...
        (firingrate.(name).contra.evoked.average - firingrate.(name).contra.baseline.average))./...
        ((firingrate.(name).ipsi.evoked.average - firingrate.(name).ipsi.baseline.average)+...
        (firingrate.(name).contra.evoked.average - firingrate.(name).contra.baseline.average));
    
end



%% Display average firing rate information
%This cell displays the coordinates of your penetration. It also displays
%the ipsi and contra evoked and basline firingrate as well as the OD score.

display(coordinates) %displays the structure coordinates
display(['Firing rate calculated from stimulus onset to ' num2str(window.firing_window) ' msecs'])
for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron
    
    display(name) %displays the neuron's name
    display(['contra evoked: ' num2str(firingrate.(name).contra.evoked.average)]) %displays the ipsi evoked FR
    display(['contra baseline: ' num2str(firingrate.(name).contra.baseline.average)]) %displays the ipsi baseline FR
    display(['contra significant? ' firingrate.(name).contra.sig_diff])
    display(['contra evoked latency (median): ' num2str(firingrate.(name).contra.evoked.median_latency)]);
    
    display(['ipsi evoked: ' num2str(firingrate.(name).ipsi.evoked.average)]) %displays the contra evoked FR
    display(['ipsi baseline: ' num2str(firingrate.(name).ipsi.baseline.average)]) %displays the contra baseline FR
    display(['ipsi significant? ' firingrate.(name).ipsi.sig_diff])
    display(['ipsi evoked latency (median): ' num2str(firingrate.(name).ipsi.evoked.median_latency)]);
    
    display(['OD score: ' num2str(firingrate.(name).ODscore)]) %displays the OD score
    display(' ') %creates blank line
    display(' ')
    
end



%% Making the structure used to construct time histograms - PSTH
%This cell makes the strucutres with spike timing information to 
%make the PSTHs. While the code looks quite complicated (because it is a 
%lot of indexing through various structures, the idea is simple)

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron
    
    %This section runs two for loops (one for ipsi and one for contra)
    %to calculate the time of spikes relative to the start of the stimulus.
    %For each trial, it takes the time at which the stimulus occured
    %(timestamps.events.ipsi or contra), and looks for all the timestamps
    %of spikes that occured in the neuron of interest between the prestimulus time
    %and poststimulus time. This is then subtracted by the time of the
    %visual stimulus, thereby setting the onset of the stimulus as time
    %zero.
    for j = 1:numberoftrials
        psth.(name).ipsi{j} = timestamps.(name)(timestamps.(name) > (timestamps.events.ipsi(j)-window.prestim-(window.window/2))...
            & timestamps.(name) < (timestamps.events.ipsi(j)+window.stimulus+window.prestim)+(window.window/2))...
            - timestamps.events.ipsi(j);
    end
    
    for j = 1:numberoftrials
        psth.(name).contra{j} = timestamps.(name)(timestamps.(name) > (timestamps.events.contra(j)-window.prestim-(window.window/2))...
            & timestamps.(name) < (timestamps.events.contra(j)+window.stimulus+window.prestim)+(window.window/2))...
            - timestamps.events.contra(j);
    end
end



%% Making the structure used to construct time histograms - spike_count
    %FOR SMOOTH PSTH
    %This nested set of loops calculates the number of spikes that occur in
    %a bin of length window.window centered around the current value of
    %centers. It bins across all trials, and therefore must be divided by
    %the number of trials to get an average spike count. Using the psth
    %structure created above, it finds the number spikes that occured 
    %around the current center with a bin of length window.window.
    
for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron
    centers = -window.prestim:window.step:(window.stimulus+window.prestim); %creates a matrix for the timing of the center of bins
    for jj = 1:length(centers)
        spike_count.(name).ipsi(jj) = 0;
        spike_count.(name).contra(jj) = 0;
        
        for ii = 1:length(psth.(name).ipsi)
            spike_count.(name).ipsi(jj) = spike_count.(name).ipsi(jj) + ((length(find(psth.(name).ipsi{ii} > (centers(jj)-(window.window/2))...
                & psth.(name).ipsi{ii} < (centers(jj) + (window.window/2)))))./((window.window/1000)*numberoftrials));
        end
        
        for ii = 1:length(psth.(name).contra)
            spike_count.(name).contra(jj) = spike_count.(name).contra(jj) + ((length(find(psth.(name).contra{ii} > (centers(jj)-(window.window/2))...
                & psth.(name).contra{ii} < (centers(jj) + (window.window/2)))))./((window.window/1000)*numberoftrials));
        end  
    end
end



%% Making smooth peristimulus time histograms
%This cell makes the smooth peristimulus time histograms for each eye. It
%will show a curve corresponding to the firing rate from window.prestim
%before the stimulus to window.prestim msecs after the stimulus. Each point
%in the curve is the number of spikes in the bin window.window. The time
%between each window.window is window.step.

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron    
    
    %for separate histograms
    figure()
    subplot(2,1,2), hold on
    title(['Smooth histogram with ' num2str(window.bin) 'ms bins of ' name ':    ' 'Ipsi'])
    db_preparepsth(window,spike_count,name) %prepares the axes, draws lines for start and stop of stimulus
    plot(centers,spike_count.(name).ipsi,'-r') %plots PSTH
    hold off
    
    subplot(2,1,1), hold on
    title(['Smooth histogram with ' num2str(window.bin) 'ms bins of ' name ':    ' 'Contra'])
    db_preparepsth(window,spike_count,name)
    plot(centers,spike_count.(name).contra,'-b')
    hold off
    saveas(figure(gcf),['D:/' user '/' folder_name  '_' date_of_exp '/Smooth_PSTH_' name '_' date_of_exp],'fig')
end



%% Making the structure used to construct time histograms - spike_count_bin
    %FOR BIN PSTH
    %This nested set of loops calculates the number of spikes that occur in
    %a bin of length window.bin from current centers_bin to the next one. 
    %It bins across all trials, and therefore must be divided by
    %the number of trials to get an average spike count.
    
for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron    
    
    centers_bin = -window.prestim:window.bin:(window.stimulus+window.prestim);
    for jj = 1:length(centers_bin)-1
        spike_count_bin.(name).ipsi(jj) = 0;
        spike_count_bin.(name).contra(jj) = 0;

        for ii = 1:length(psth.(name).ipsi)
            spike_count_bin.(name).ipsi(jj) = spike_count_bin.(name).ipsi(jj) +...
                ((length(find(psth.(name).ipsi{ii} > centers_bin(jj) & psth.(name).ipsi{ii} < (centers_bin(jj+1)))))./((window.bin/1000)*numberoftrials));
        end
        
        for ii = 1:length(psth.(name).contra)
            spike_count_bin.(name).contra(jj) = spike_count_bin.(name).contra(jj) +...
                ((length(find(psth.(name).contra{ii} > centers_bin(jj) & psth.(name).contra{ii} < (centers_bin(jj+1)))))./((window.bin/1000)*numberoftrials));
        end
    end
    
end




%% Make bin histogram PSTH

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron 
    
    %allows for the histogram to be centered in the middle of each
    %window.bin
    centered_for_hist = -window.prestim+window.bin/2:window.bin:(window.stimulus+window.prestim);
    
    %for separate histograms
    figure()
    subplot(2,1,2), hold on
    title(['Histogram with ' num2str(window.bin) 'ms bins of ' name ':    ' 'Ipsi'])
    db_preparepsth(window,spike_count_bin,name)
    bar(centered_for_hist,spike_count_bin.(name).ipsi, 'r') %plots the histogram
    hold off
    
    subplot(2,1,1), hold on
    title(['Histogram with ' num2str(window.bin) 'ms bins of ' name ':    ' 'Contra'])
    db_preparepsth(window,spike_count_bin,name)
    bar(centered_for_hist,spike_count_bin.(name).contra, 'b')
    hold off
    saveas(figure(gcf),['D:/' user '/' folder_name  '_' date_of_exp '/Bin_PSTH_' name '_' date_of_exp],'fig')
end

%% Make rasters

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)];
    
    %
    figure()
    
    for k = 1:length(eye_open)
        subplot(2,1,k), hold on
        title(['Rasters for ' name ':  ' eye_open{k}])
        ylim([0 numberoftrials]);
        ylabel('Trial')
        xlim([-window.prestim window.stimulus+window.prestim]);
        xlabel('Time (msec)')
        
        %draw lines for when stimulus turned on and off
        line([0 0],get(gca, 'YLim'), 'Color', [.5 .5 .5], 'LineWidth', 2, 'LineStyle', '-')
        line([window.stimulus window.stimulus], get(gca, 'YLim'), 'Color', [.5 .5 .5], 'LineWidth', 2, 'LineStyle', '--')
        
        %draws line for each spike in raster
        for j = 1:numberoftrials
            for ii = 1:length(firingrate.(name).(eye_open{k}).evoked.for_rasters{j})
                line([firingrate.(name).(eye_open{k}).evoked.for_rasters{j}(ii)...
                    firingrate.(name).(eye_open{k}).evoked.for_rasters{j}(ii)],...
                    [j-1 j], 'LineWidth', 1, 'Color', [0 0 0])
            end
        end
    end
    
    saveas(figure(gcf),['D:/' user '/' folder_name  '_' date_of_exp '/Raster_' name '_' date_of_exp],'fig')
end

%% Calculating Firing Rate using a sliding window (window.window)
%This will find the max firing rate in the smooth PSTH and will go back
%through the dataset by trial to give you the standard deviation. The
%baseline is calculated using a same size window just before the 
%presentation of the stimulus. OD scores are calculated by the same 
%formula stated in 'Calculate Firing Rate.'

%A paired t-test is used to assess whether the contra or ipsi evoked
%response is signficantly above baseline.

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)];
    
    %loop for either ipsi or contra eye
    for k = 1:length(eye_open)
        %tells you the size of the window (msecs)
        max_fr.(name).(eye_open{k}).firing_rate = window.window;
        
        %finds the maximum firing rate during the stimulus
        max_fr.(name).(eye_open{k}).firing_rate = max(spike_count.(name).(eye_open{k})(centers >= 0 & centers <= window.stimulus));
        
        %finds the all the times this firing rate is reached
        max_fr.(name).(eye_open{k}).time_max =...
            centers(spike_count.(name).(eye_open{k}) == max(spike_count.(name).(eye_open{k})(centers >= 0 & centers <= window.stimulus)));
        %gets rid of any time this firing rate is reached before the
        %stimulus
        max_fr.(name).(eye_open{k}).time_max = max_fr.(name).(eye_open{k}).time_max(max_fr.(name).(eye_open{k}).time_max > 0);
        %finds the earliest time after the stimulus when the FR is reached
        max_fr.(name).(eye_open{k}).time_max = min(max_fr.(name).(eye_open{k}).time_max);
        
        %checks to see if time to max FR is shorter than half the
        %window.window size. If so, time is set to window.window/2.
        %Prevents prestim spikes from being counted.
        if max_fr.(name).(eye_open{k}).time_max >= window.window/2
            max_fr.(name).(eye_open{k}).time_used = max_fr.(name).(eye_open{k}).time_max;
        else
            max_fr.(name).(eye_open{k}).time_used = window.window/2;
        end
        
        max_fr.(name).(eye_open{k}).per_trial = zeros(1,numberoftrials);
        max_fr.(name).(eye_open{k}).spontaneous_trial = zeros(1,numberoftrials);
        for j = 1:numberoftrials
            max_fr.(name).(eye_open{k}).per_trial(j) =...
                sum(psth.(name).(eye_open{k}){j} > max_fr.(name).(eye_open{k}).time_used - window.window/2 &...
                psth.(name).(eye_open{k}){j} < max_fr.(name).(eye_open{k}).time_used + window.window/2);
            
            max_fr.(name).(eye_open{k}).spontaneous_trial(j) = ...
                sum(psth.(name).(eye_open{k}){j} > 0 - window.window &...
                psth.(name).(eye_open{k}){j} < 0);
        end
        
        %Calculates sd of maximal response
        max_fr.(name).(eye_open{k}).sd = std(max_fr.(name).(eye_open{k}).per_trial);
        
        %Baseline firing rate and sd
        max_fr.(name).(eye_open{k}).baseline_FR = mean(max_fr.(name).(eye_open{k}).spontaneous_trial);
        max_fr.(name).(eye_open{k}).baseline_SD = std(max_fr.(name).(eye_open{k}).spontaneous_trial);
        
        %Is maxiaml response signficantly different than baseline?
        [max_fr.(name).(eye_open{k}).paired_test.h, max_fr.(name).(eye_open{k}).paired_test.p, ...
            max_fr.(name).(eye_open{k}).paired_test.CI] = ...
            ttest(max_fr.(name).(eye_open{k}).per_trial,max_fr.(name).(eye_open{k}).spontaneous_trial);
        
        if max_fr.(name).(eye_open{k}).paired_test.h == 1
            max_fr.(name).(eye_open{k}).diff_from_base = 'yes';
        elseif max_fr.(name).(eye_open{k}).paired_test.h == 0
            max_fr.(name).(eye_open{k}).diff_from_base = 'no';
        end
        
    end
    
    %calculates OD score
    % ODscore = [(ipsi evoked - ipsi baseline) - (contra evoked - contra baseline)] / [(ipsi evoked - ipsi baseline) + (contra evoked - contra baseline)]
    % Score of -1 is all contra, score of 1 is all ipsi
    
    max_fr.(name).OD = ((max_fr.(name).ipsi.firing_rate - max_fr.(name).ipsi.baseline_FR) -...
        (max_fr.(name).contra.firing_rate - max_fr.(name).contra.baseline_FR))./...
        ((max_fr.(name).ipsi.firing_rate - max_fr.(name).ipsi.baseline_FR) +...
        (max_fr.(name).contra.firing_rate - max_fr.(name).contra.baseline_FR));
end
    
    
%% Display Results from sliding window firing rate calculations                
        
display(' ')
display(['Maximal firing rate calculated with a window of ' num2str(window.window) ' msecs'])

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron
    
    display(name) %displays the neuron's name
    display(['contra evoked: ' num2str(max_fr.(name).contra.firing_rate)]) %displays the contra evoked FR
    display(['contra baseline: ' num2str(max_fr.(name).contra.baseline_FR)]) %displays the contra baseline FR
    display(['contra significant? ' max_fr.(name).contra.diff_from_base]) %says if contra FR is signifcantly above baseline
    
    display(['ipsi evoked: ' num2str(max_fr.(name).ipsi.firing_rate)])
     display(['ipsi baseline: ' num2str(max_fr.(name).ipsi.baseline_FR)]) 
    display(['ipsi significant? ' max_fr.(name).ipsi.diff_from_base])
    
    display(['OD score: ' num2str(max_fr.(name).OD)]) %displays the OD score
    display(' ') %creates blank line
    display(' ')
    
end

%% Save workspace
%Will save all the workspace variables in case you want to look at it
%again.

save(['D:/' user '/' folder_name '_' date_of_exp '/workspace' '_' folder_name '_' date_of_exp])

    
    
    








    
    
    