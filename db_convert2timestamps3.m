function [ timestamps ] = db_convert2timestamps3(data, header, varargin)
%db_convert2timestamps3(data, header, {time_between_opto}) Converts a 'data' 
%matrix of {0,1} spike times/events to a series of timestamps. Uses information 
%from the cell 'header' to determine the number of neurons, their channel, and
%the type of events (visual and/or optogenetic stimulation). Also asks which 
%visual events are contra or ipsi. Can take an additional argument about the 
%timing between optostimulation bouts. Will be 1000ms by default.
%
%Written by DM Brady, November 2013.


%Searches for a column in header called EVT1 (meaning visual events)
if sum(strcmpi(header,'EVT1')) == 1
  %finds the column in data that has the timing of visual events
  timestamps.visual_events.col_num = find(strcmpi(header,'EVT1') == 1);
  %finds timestamps for all visual events
  timestamps.visual_events.all = find(data(:,timestamps.visual_events.col_num) > 0);

  %Asks for the trial #'s for the contra eye
  contra =  input(['What trials were for the contra eye?\n There are a total of ' ...
    num2str(length(timestamps.visual_events.all)) ...
    ' events \n ex: [1 2 3] [1:4 9:12] [1:2:20]??   ']);

  %Asks for the trial #'s for the ipsi eye
  ipsi =  input(['What trials were for the ipsi eye?\n There are a total of ' ...
    num2str(length(timestamps.visual_events.all)) ...
    ' events \n ex: [1 2 3] [1:4 9:12] [1:2:20]??   ']);

  timestamps.visual_events.contra = timestamps.visual_events.all(contra); %event times for contra eye
  timestamps.visual_events.ipsi = timestamps.visual_events.all(ipsi); %event times for the ipsi eye
end

%Searches for a column in header called EVT17 (meaning optogenetic stimulation)
if sum(strcmpi(header,'EVT17')) == 1
  %finds the column in data that has the timing of optogenetic events
  timestamps.opto_events.col_num = find(strcmpi(header,'EVT17') == 1);
  %finds the timestamps for all optogenetic events
  timestamps.opto_events.all = find(data(:,timestamps.opto_events.col_num) > 0);

  %finds the timestamps for first optogenetic stimulation in train
  %looks for an argument for the time between bouts, if it is not specified,
  %it is set to 1000 ms.
  switch nargin
    case 2
      time_between_bouts = 1000;
    otherwise
      time_between_bouts = varargin{1};
  end
  
  %takes the first opto event and every opto event after a time_between_bouts
  timestamps.opto_events.first = [timestamps.opto_events.all(1);...
    timestamps.opto_events.all(circshift(diff(timestamps.opto_events.all)>time_between_bouts,1))];
end

%Searches for number of neurons
timestamps.neuron_sum.num = sum(~cellfun('isempty',strfind(header,'SPK')));
timestamps.neuron_sum.indices = find(~cellfun('isempty',strfind(header,'SPK')));
for i = 1:timestamps.neuron_sum.num
  timestamps.(['neuron_' num2str(i)]).timestamps = find(data(:,timestamps.neuron_sum.indices(i)) > 0);
  timestamps.(['neuron_' num2str(i)]).channel = str2double(header{timestamps.neuron_sum.indices(i)}(4:5));
end



end
