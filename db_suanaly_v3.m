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


parameters.user = input('What is your name?  ', 's'); %asks for your name to make a folder

parameters.folder_name = input('What is the name of your experiment?  ', 's'); %asks for the name of your experiment
parameters.date_of_exp = input('When did you do this experiment? \n(example: 06Jul2012)  ', 's'); %asks for the date of your experiment
parameters.penetration_number = ['pen_' input('What penetration is this? ','s')];

%coordinates will be a structure with three different fields corresponding
%to the A/P, M/L, and depth of your recording. Type coordinate.(field of
%interest) to look at that number.
parameters.coordinates.AP = input('What are your anterior/posterior coordinates?  '); %asks for the anterior/posterior coordinates
parameters.coordinates.ML = input('What are your medial/lateral coordinates?  '); %asks for the medial/lateral coordinates
parameters.coordinates.depth = input('What is the depth of your probe?  '); %asks for the depth of your probe

parameters.folder_save = ['D:/' parameters.user '/' parameters.folder_name '/' parameters.date_of_exp];
mkdir([parameters.folder_save '/' parameters.penetration_number]) %the actual function
%that will create the folder to store all of your information. Currently,
%it will save in the D:/ drive. Ex:
%  D:/(user)/(folder name)_(date of experiment)
%  D:/Daniel/LynxKO/05Jul2012/pen1
if exist([parameters.folder_save '/good_cells'],'dir') == 0
    mkdir([parameters.folder_save '/good_cells'])
end

if exist([parameters.folder_save '/good_cells/good_cells.csv'],'file') == 0
    col_header = {'Contra_evoked','Contra_base','Contra_delta',...
      'Ipsi_evoked','Ipsi_base','Ipsi_delta','ODI','ODS','SU/MU'};
    [nrow, ncol] = size(col_header);
    fid = fopen([parameters.folder_save '/good_cells/good_cells.csv'],'w');
    for i = 1:ncol
        if i < ncol
            fprintf(fid,'%s,',col_header{:,i});
        else
            fprintf(fid,'%s\n',col_header{:,i});
        end
    end
    fclose(fid);
end



%% Convert single unit data from plexon into a structure of 'timestamps'
%This cell will calculate the number of visual stimuli used and convert
%data from a sequence of 0's and 1's to a structure with the times at which
%events occured (whether they are spikes or the presentation fo a visual
%stimulus).

[timestamps, parameters.numberoftrials] = ...
 db_convert2timestamps2(data); %Calls the function db_convert2timestamps2 which transforms
%the variable 'data' from your recording to the times when spikes or visual
%stimuli occured.



%% Create window structure for calculating firing rates and making peristimulus histograms
%The structure window will have values for the stimulus length, pre and
%post stim period for making histograms, and the size of bins for the
%smooth and binned PSTH

%Stimulus length and pre/post stimulus period
parameters.window.stimulus = 5000; %number of milliseconds the stimulus is on
parameters.window.prestim = 2000; %number of milliseconds before stimulus to calculate baseline (also used to plot post stimulus spikes)

%Smooth PSTH 
parameters.window.window = 200; %length (in msec) of bin for psth
parameters.window.step = 2; %time difference between parameters.centers of windows

%Maximal Firing Rate Window
parameters.window.fr = 500;

%Bin PSTH
parameters.window.bin = 200; %length (in msec) of bin for histogram


parameters.eye_open = {'contra' 'ipsi'};

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
    for j = 1:parameters.numberoftrials
        psth.(name).ipsi{j} = timestamps.(name)(timestamps.(name) > (timestamps.events.ipsi(j)-parameters.window.prestim-(parameters.window.window/2))...
            & timestamps.(name) < (timestamps.events.ipsi(j)+parameters.window.stimulus+parameters.window.prestim)+(parameters.window.window/2))...
            - timestamps.events.ipsi(j);
    end
    
    for j = 1:parameters.numberoftrials
        psth.(name).contra{j} = timestamps.(name)(timestamps.(name) > (timestamps.events.contra(j)-parameters.window.prestim-(parameters.window.window/2))...
            & timestamps.(name) < (timestamps.events.contra(j)+parameters.window.stimulus+parameters.window.prestim)+(parameters.window.window/2))...
            - timestamps.events.contra(j);
    end
end



%% Making the structure used to construct time histograms - spike_count
    %FOR SMOOTH PSTH
    %This nested set of loops calculates the number of spikes that occur in
    %a bin of length window.window centered around the current value of
    %parameters.centers. It bins across all trials, and therefore must be divided by
    %the number of trials to get an average spike count. Using the psth
    %structure created above, it finds the number spikes that occured 
    %around the current center with a bin of length window.window.
    
for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron
    parameters.centers = -parameters.window.prestim:parameters.window.step: ...
     (parameters.window.stimulus+parameters.window.prestim); %creates a matrix for the timing of the center of bins
    for jj = 1:length(parameters.centers)
        spike_count.(name).ipsi(jj) = 0;
        spike_count.(name).contra(jj) = 0;
        
        for ii = 1:length(psth.(name).ipsi)
            spike_count.(name).ipsi(jj) = spike_count.(name).ipsi(jj) + ((length(find(psth.(name).ipsi{ii} > (parameters.centers(jj)-(parameters.window.window/2))...
                & psth.(name).ipsi{ii} < (parameters.centers(jj) + (parameters.window.window/2)))))./((parameters.window.window/1000)*parameters.numberoftrials));
        end
        
        for ii = 1:length(psth.(name).contra)
            spike_count.(name).contra(jj) = spike_count.(name).contra(jj) + ((length(find(psth.(name).contra{ii} > (parameters.centers(jj)-(parameters.window.window/2))...
                & psth.(name).contra{ii} < (parameters.centers(jj) + (parameters.window.window/2)))))./((parameters.window.window/1000)*parameters.numberoftrials));
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
    figure(i), title(name)
    subplot(4,2,3), hold on
    title('Ipsi PSTH')
    db_preparepsth(parameters.window,spike_count,name) %prepares the axes, draws lines for start and stop of stimulus
    plot(parameters.centers,spike_count.(name).ipsi,'-r') %plots PSTH
    hold off
    
    subplot(4,2,1), hold on
    title('Contra PSTH')
    db_preparepsth(parameters.window,spike_count,name)
    plot(parameters.centers,spike_count.(name).contra,'-b')
    hold off
end



%% Making the structure used to construct time histograms - spike_count_bin
    %FOR BIN PSTH
    %This nested set of loops calculates the number of spikes that occur in
    %a bin of length window.bin from current parameters.centers_bin to the next one. 
    %It bins across all trials, and therefore must be divided by
    %the number of trials to get an average spike count.
    
for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron    
    
    parameters.centers_bin = -parameters.window.prestim:parameters.window.bin:(parameters.window.stimulus+parameters.window.prestim);
    for jj = 1:length(parameters.centers_bin)-1
        spike_count_bin.(name).ipsi(jj) = 0;
        spike_count_bin.(name).contra(jj) = 0;

        for ii = 1:length(psth.(name).ipsi)
            spike_count_bin.(name).ipsi(jj) = spike_count_bin.(name).ipsi(jj) +...
                ((length(find(psth.(name).ipsi{ii} > parameters.centers_bin(jj) & psth.(name).ipsi{ii} < (parameters.centers_bin(jj+1)))))./((parameters.window.bin/1000)*parameters.numberoftrials));
        end
        
        for ii = 1:length(psth.(name).contra)
            spike_count_bin.(name).contra(jj) = spike_count_bin.(name).contra(jj) +...
                ((length(find(psth.(name).contra{ii} > parameters.centers_bin(jj) & psth.(name).contra{ii} < (parameters.centers_bin(jj+1)))))./((parameters.window.bin/1000)*parameters.numberoftrials));
        end
    end
    
end




%% Make bin histogram PSTH

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron 
    
    %allows for the histogram to be centered in the middle of each
    %window.bin
    parameters.centered_for_hist = -parameters.window.prestim+parameters.window.bin/2:parameters.window.bin:(parameters.window.stimulus+parameters.window.prestim);
    
    %for separate histograms
    figure(i)
    subplot(4,2,4), hold on
    title(['Ipsi PSTH'])
    db_preparepsth(parameters.window,spike_count_bin,name)
    bar(parameters.centered_for_hist,spike_count_bin.(name).ipsi, 'r') %plots the histogram
    hold off
    
    subplot(4,2,2), hold on
    title(['Contra PSTH'])
    db_preparepsth(parameters.window,spike_count_bin,name)
    bar(parameters.centered_for_hist,spike_count_bin.(name).contra, 'b')
    hold off
    
end

%% Make rasters
for i = 1:size(data,2)-1 %to repeat the following loop for the number of neurons recorded in this penetration
    
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron
    
    %loop for either ipsi or contra eye
    for k = 1:length(parameters.eye_open)
    
        for j = 1:parameters.numberoftrials
            firingrate.(name).(parameters.eye_open{k}).evoked.for_rasters{j} = find(data(timestamps.events.(parameters.eye_open{k})(j)-parameters.window.prestim...
                :timestamps.events.(parameters.eye_open{k})(j)...
                +parameters.window.stimulus+parameters.window.prestim,i))-parameters.window.prestim;
        end
    end
end
                
                
for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)];
    
    %
    figure(i)
    
    for k = 1:length(parameters.eye_open)
        subplot(4,2,2*k+3), hold on
        title([parameters.eye_open{k} ' Rasters'])
        ylim([0 parameters.numberoftrials]);
        ylabel('Trial')
        xlim([-parameters.window.prestim parameters.window.stimulus+parameters.window.prestim]);
        xlabel('Time (msec)')
        
        %draw lines for when stimulus turned on and off
        line([0 0],get(gca, 'YLim'), 'Color', [.5 .5 .5], 'LineWidth', 2, 'LineStyle', '-')
        line([parameters.window.stimulus parameters.window.stimulus], get(gca, 'YLim'), 'Color', [.5 .5 .5], 'LineWidth', 2, 'LineStyle', '--')
        
        %draws line for each spike in raster
        for j = 1:parameters.numberoftrials
            for ii = 1:length(firingrate.(name).(parameters.eye_open{k}).evoked.for_rasters{j})
                line([firingrate.(name).(parameters.eye_open{k}).evoked.for_rasters{j}(ii)...
                    firingrate.(name).(parameters.eye_open{k}).evoked.for_rasters{j}(ii)],...
                    [j-1 j], 'LineWidth', 1, 'Color', [0 0 0])
            end
        end
    end
end



%% Calculating Firing Rate using a sliding window (window.window)
%This will find the max firing rate in the smooth PSTH and will go back
%through the dataset by trial to give you the standard deviation. The
%baseline is calculated using a same size window just before the 
%presentation of the stimulus. OD scores are calculated by the same 
%formula stated in 'Calculate Firing Rate.'


for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)]; %name of the 'ith' neuron
    parameters.centers = -parameters.window.prestim:parameters.window.step:(parameters.window.stimulus+parameters.window.prestim); %creates a matrix for the timing of the center of bins
    for jj = 1:length(parameters.centers)
        spike_count_fr.(name).ipsi(jj) = 0;
        spike_count_fr.(name).contra(jj) = 0;
        
        for ii = 1:length(psth.(name).ipsi)
            spike_count_fr.(name).ipsi(jj) = spike_count_fr.(name).ipsi(jj) + ((length(find(psth.(name).ipsi{ii} > (parameters.centers(jj)-(parameters.window.fr/2))...
                & psth.(name).ipsi{ii} < (parameters.centers(jj) + (parameters.window.fr/2)))))./((parameters.window.fr/1000)*parameters.numberoftrials));
        end
        
        for ii = 1:length(psth.(name).contra)
            spike_count_fr.(name).contra(jj) = spike_count_fr.(name).contra(jj) + ((length(find(psth.(name).contra{ii} > (parameters.centers(jj)-(parameters.window.fr/2))...
                & psth.(name).contra{ii} < (parameters.centers(jj) + (parameters.window.fr/2)))))./((parameters.window.fr/1000)*parameters.numberoftrials));
        end  
    end
end

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)];
    
    %loop for either ipsi or contra eye
    for k = 1:length(parameters.eye_open)
        %tells you the size of the window (msecs)
        max_fr.(name).(parameters.eye_open{k}).firing_rate = parameters.window.fr;
        
        %finds the maximum firing rate during the stimulus
        max_fr.(name).(parameters.eye_open{k}).firing_rate = max(spike_count_fr.(name).(parameters.eye_open{k})(parameters.centers >= 0 & parameters.centers <= parameters.window.stimulus));
        
        %finds the all the times this firing rate is reached
        max_fr.(name).(parameters.eye_open{k}).time_max =...
            parameters.centers(spike_count_fr.(name).(parameters.eye_open{k}) == max(spike_count_fr.(name).(parameters.eye_open{k})(parameters.centers >= 0 & parameters.centers <= parameters.window.stimulus)));
        %gets rid of any time this firing rate is reached before the
        %stimulus
        max_fr.(name).(parameters.eye_open{k}).time_max = max_fr.(name).(parameters.eye_open{k}).time_max(max_fr.(name).(parameters.eye_open{k}).time_max > 0);
        %finds the earliest time after the stimulus when the FR is reached
        max_fr.(name).(parameters.eye_open{k}).time_max = min(max_fr.(name).(parameters.eye_open{k}).time_max);
        
        %checks to see if time to max FR is shorter than half the
        %window.window size. If so, time is set to window.fr/2.
        %Prevents prestim spikes from being counted.
        if max_fr.(name).(parameters.eye_open{k}).time_max >= parameters.window.fr/2
            max_fr.(name).(parameters.eye_open{k}).time_used = max_fr.(name).(parameters.eye_open{k}).time_max;
        else
            max_fr.(name).(parameters.eye_open{k}).time_used = parameters.window.fr/2;
        end
        
        max_fr.(name).(parameters.eye_open{k}).per_trial = zeros(1,parameters.numberoftrials);
        %max_fr.(name).(parameters.eye_open{k}).spontaneous_trial = zeros(1,parameters.numberoftrials);
        
        for j = 1:parameters.numberoftrials
            max_fr.(name).(parameters.eye_open{k}).per_trial(j) =...
                sum(psth.(name).(parameters.eye_open{k}){j} > max_fr.(name).(parameters.eye_open{k}).time_used - parameters.window.fr/2 &...
                psth.(name).(parameters.eye_open{k}){j} < max_fr.(name).(parameters.eye_open{k}).time_used + parameters.window.fr/2);
            
            %max_fr.(name).(parameters.eye_open{k}).spontaneous_trial(j) = ...
             %   sum(psth.(name).(parameters.eye_open{k}){j} > 0 - parameters.window.fr &...
              %  psth.(name).(parameters.eye_open{k}){j} < 0);
        end
        
        % For calculating spontaneous response (bootstrap with same window size)
        boot_num = 10000 * length(timestamps.events.(parameters.eye_open{k}));
        trial_selection = randi(length(timestamps.events.(parameters.eye_open{k})),boot_num,1);
        time_selection = -1*round((parameters.window.prestim-parameters.window.fr)*rand(boot_num,1)+parameters.window.fr);
        for j = 1:boot_num
            max_fr.(name).(parameters.eye_open{k}).spontaneous.boot(j) = sum(psth.(name).(parameters.eye_open{k}){trial_selection(j)} > ...
                time_selection(j) & psth.(name).(parameters.eye_open{k}){trial_selection(j)} < time_selection(j)+parameters.window.fr);
        end
        
        max_fr.(name).(parameters.eye_open{k}).spontaneous.boot =...
         reshape(max_fr.(name).(parameters.eye_open{k}).spontaneous.boot, boot_num/length(timestamps.events.(parameters.eye_open{k})),...
         length(timestamps.events.(parameters.eye_open{k})));
        max_fr.(name).(parameters.eye_open{k}).spontaneous.mean = ...
         median(mean(max_fr.(name).(parameters.eye_open{k}).spontaneous.boot,2))./(parameters.window.fr/1000);
        max_fr.(name).(parameters.eye_open{k}).spontaneous.CI = [prctile(mean(max_fr.(name).(parameters.eye_open{k}).spontaneous.boot,2),97.5) ...
         prctile(mean(max_fr.(name).(parameters.eye_open{k}).spontaneous.boot,2),2.5)]./(parameters.window.fr/1000);
        
        %Calculates mean, sd of maximal response
        %max_fr.(name).(parameters.eye_open{k}).mean = mean(max_fr.(name).(parameters.eye_open{k}).per_trial);
        %max_fr.(name).(parameters.eye_open{k}).sd = std(max_fr.(name).(parameters.eye_open{k}).per_trial);
        
        %Baseline firing rate and sd
        %max_fr.(name).(parameters.eye_open{k}).baseline_FR = mean(max_fr.(name).(parameters.eye_open{k}).spontaneous_trial);
        %max_fr.(name).(parameters.eye_open{k}).baseline_SD = std(max_fr.(name).(parameters.eye_open{k}).spontaneous_trial);
        
        %Is maximal response signficantly different than baseline?
        %[max_fr.(name).(parameters.eye_open{k}).paired_test.h, max_fr.(name).(parameters.eye_open{k}).paired_test.p, ...
         %   max_fr.(name).(parameters.eye_open{k}).paired_test.CI] = ...
          %  ttest(max_fr.(name).(parameters.eye_open{k}).per_trial,max_fr.(name).(parameters.eye_open{k}).spontaneous_trial);
        
        if max_fr.(name).(parameters.eye_open{k}).firing_rate > max_fr.(name).(parameters.eye_open{k}).spontaneous.CI(1)
            max_fr.(name).(parameters.eye_open{k}).diff_from_base = 'yes';
        elseif max_fr.(name).(parameters.eye_open{k}).firing_rate < max_fr.(name).(parameters.eye_open{k}).spontaneous.CI(2)
            max_fr.(name).(parameters.eye_open{k}).diff_from_base = 'yes';
        else
            max_fr.(name).(parameters.eye_open{k}).diff_from_base = 'no';
        end
        
      %Calculates normalized firing rate (stimulus - blank)/(stimulus + blank)
      max_fr.(name).(parameters.eye_open{k}).normal_fr = (max_fr.(name).(parameters.eye_open{k}).firing_rate - ...
         max_fr.(name).(parameters.eye_open{k}).spontaneous.mean) ./ ...
         (max_fr.(name).(parameters.eye_open{k}).firing_rate + ...
          max_fr.(name).(parameters.eye_open{k}).spontaneous.mean);

    end
    
    %calculates OD score
    % ODscore = [(ipsi evoked - ipsi baseline) - (contra evoked - contra baseline)] / [(ipsi evoked - ipsi baseline) + (contra evoked - contra baseline)]
    % Score of -1 is all contra, score of 1 is all ipsi
    
    max_fr.(name).ODI = ((max_fr.(name).ipsi.firing_rate - max_fr.(name).ipsi.spontaneous.mean) -...
        (max_fr.(name).contra.firing_rate - max_fr.(name).contra.spontaneous.mean))./...
        ((max_fr.(name).ipsi.firing_rate - max_fr.(name).ipsi.spontaneous.mean) +...
        (max_fr.(name).contra.firing_rate - max_fr.(name).contra.spontaneous.mean));

end

%Calculate ODS from ODI

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)];
    
    if max_fr.(name).ODI <= -.5
        max_fr.(name).ODS = 1;
    elseif max_fr.(name).ODI > -.5 && max_fr.(name).ODI <= -.3
        max_fr.(name).ODS = 2;
    elseif max_fr.(name).ODI > -.3 && max_fr.(name).ODI <= -.1
        max_fr.(name).ODS = 3;
    elseif max_fr.(name).ODI > -.1 && max_fr.(name).ODI <= .1
        max_fr.(name).ODS = 4;
    elseif max_fr.(name).ODI > .1 && max_fr.(name).ODI <= .3
        max_fr.(name).ODS = 5;
    elseif max_fr.(name).ODI > .3 && max_fr.(name).ODI <= .5
        max_fr.(name).ODS = 6;
    elseif max_fr.(name).ODI > .5
        max_fr.(name).ODS = 7;
    end
end
    
    
%% Display Results from sliding window firing rate calculations

for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)];
    
    %
    figure(i)
    hold on
    annotation('textbox',[.55 .25 .3 .2],'FitBoxtoText','on','String',...
        {['contra evoked: ' num2str(max_fr.(name).contra.firing_rate)];...
        ['contra baseline: ' num2str(max_fr.(name).contra.spontaneous.mean)];...
        ['contra normalized change in firing rate: '  num2str(max_fr.(name).contra.normal_fr)];...
        ['contra significant? ' max_fr.(name).contra.diff_from_base];...
        ['ipsi evoked: ' num2str(max_fr.(name).ipsi.firing_rate)];...
        ['ipsi baseline: ' num2str(max_fr.(name).ipsi.spontaneous.mean)];...
        ['ipsi normalized change in firing rate: '  num2str(max_fr.(name).ipsi.normal_fr)];...
        ['ipsi significant? ' max_fr.(name).ipsi.diff_from_base];...
        ['ODI score: ' num2str(max_fr.(name).ODI)];...
        ['ODS score: ' num2str(max_fr.(name).ODS)]});
    hold off
    saveas(figure(gcf),[parameters.folder_save '/' parameters.penetration_number '/Summary_' name ],'fig')
    saveas(figure(gcf),[parameters.folder_save '/' parameters.penetration_number '/Summary_' name ],'pdf')
    saveas(figure(gcf),[parameters.folder_save '/' parameters.penetration_number '/Summary_' name ],'epsc')
end             

%% Asks which units were good and saves that information separately
for i = 1:size(data,2)-1
    name = ['neuron' num2str(i)];
    
    good_neuron = input(['Is ' name ' a good neuron? (y or n) '],'s');
    %if strcmpi(max_fr.(name).contra.diff_from_base,'yes') & strcmpi(max_fr.(name).ipsi.diff_from_base,'yes')
    %    good_neuron = 'y';
    %else
    %    good_neuron = 'n';
    %end
    
    if strcmpi(good_neuron, 'y') == 1
        su_or_mu = input(['Is it single or multiunit? (su or mu) '], 's');
        
        %analysis.max_fr = max_fr.(name);
        %save([parameters.folder_save '/good_cells/' name '_' su_or_mu],'analysis')
        to_append = {max_fr.(name).contra.firing_rate, max_fr.(name).contra.spontaneous.mean, ...
          max_fr.(name).contra.normal_fr, max_fr.(name).ipsi.firing_rate, max_fr.(name).ipsi.spontaneous.mean, ...
          max_fr.(name).ipsi.normal_fr, max_fr.(name).ODI, max_fr.(name).ODS, su_or_mu};
        
        fid = fopen([parameters.folder_save '/good_cells/good_cells.csv'],'a');
        [nrow,ncol] = size(to_append);
        for j = 1:ncol
            if j < ncol
                fprintf(fid,'%f,',to_append{:,j});
            else
                fprintf(fid,'%s\n',to_append{:,j});
            end
        end
    end
end

%% Save workspace
%Will save all the workspace variables in case you want to look at it
%again.

save([parameters.folder_save '/' parameters.penetration_number '/workspace'])







    
    
    
