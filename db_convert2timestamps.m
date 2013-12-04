function [ timestamps ] = db_convert2timestamps(data,numberoftrials)
%convert2timestamps Converts large matrix of 1's and 0's from plexon to
%just when spikes and events occured
%Written by DM Brady 07/17/2012
%Updated 10/26/2012

%asks which eye was open first
which_eye = input('which eye is open first? (ipsi or contra) ','s');

%makes a cell with denoting which eye was open first.
if strcmpi(which_eye, 'ipsi') == 1
    eye_order = {'ipsi' 'contra'};
elseif strcmpi(which_eye, 'contra') == 1;
    eye_order = {'contra' 'ipsi'};
else
    display('Sorry, I did not understand you, please run the program again')
    return
end

%Find when in recording events occured (all)
timestamps.events.all = find(data(:,end) > 0);
timestamps.events.(eye_order{1}) = timestamps.events.all(1:numberoftrials); %events only for first open eye
timestamps.events.(eye_order{2}) = timestamps.events.all(numberoftrials+1:2*numberoftrials); %events only for second open eye

for i = 1:size(data,2)-1
  name = ['neuron' num2str(i)];
  timestamps.(name) = find(data(:,i) >= 1); %for loop that creates a cell with number of neurons and every time a spike occurs
  for j = 1:size(timestamps.(name),1)-1
      timestamps.([name '_isi'])(j,:) = [timestamps.(name)(j+1)-timestamps.(name)(j) timestamps.(name)(j)];
      timestamps.([name '_instan'])(j,:) = [1./timestamps.([name '_isi'])(j)*1000 timestamps.(name)(j)];
  end
end



end

