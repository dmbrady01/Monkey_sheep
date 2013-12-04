function [ timestamps, numberoftrials ] = db_convert2timestamps2(data)
%convert2timestamps Converts large matrix of 1's and 0's from plexon to
%just when spikes and events occured
%Written by DM Brady 01/30/2013

%Find when in recording events occured (all)
timestamps.events.all = find(data(:,end) > 0);

%Asks for the trial #'s for the contra eye
contra = input(['What trials were for the contra eye?\n There are a total of ' num2str(length(timestamps.events.all)) ...
    ' events \n ex: [1 2 3] [1:4 9:12] [1:2:20]??   ']);

%Asks for the trial #'s for the ipsi eye
ipsi = input(['What trials were for the ipsi eye?\n There are a total of ' num2str(length(timestamps.events.all)) ...
    ' events \n ex: [1 2 3] [1:4 9:12] [1:2:20]??   ']);

%Calculates the number of trials (if there are different numbers of trials
%put into contra and ipsi it will eliminate later trials of the larger
%condition)
if length(contra) >= length(ipsi)
    numberoftrials = length(ipsi);
elseif length(contra) <= length(ipsi)
    numberoftrials = length(contra);
end

%Find when in recording events occured (all)
timestamps.events.all = find(data(:,end) > 0);

timestamps.events.contra = timestamps.events.all(contra); %events only for contra eye
timestamps.events.ipsi = timestamps.events.all(ipsi); %events only for ipsi eye

for i = 1:size(data,2)-1
  name = ['neuron' num2str(i)];
  timestamps.(name) = find(data(:,i) >= 1); %for loop that creates a cell with number of neurons and every time a spike occurs
  for j = 1:size(timestamps.(name),1)-1
      timestamps.([name '_isi'])(j,:) = [timestamps.(name)(j+1)-timestamps.(name)(j) timestamps.(name)(j)];
      timestamps.([name '_instan'])(j,:) = [1./timestamps.([name '_isi'])(j)*1000 timestamps.(name)(j)];
  end
end



end