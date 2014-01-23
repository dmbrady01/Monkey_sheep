% db_optoanaly
%
% Used to make rasters and psth to analyze optogenetic stimulation of different cell types
%
% Written by DM Brady 12/2013

% Time between bouts of stimulation
optoanaly.parameters.btw_bouts = 150; %1000ms between opto stimulation bouts

%% Gather timestamp information of optogenetic stimulation and spiking
optoanaly.timestamps = db_convert2timestamps3(data, nexColumnNames, optoanaly.parameters.btw_bouts);

%% Setting window parameters
optoanaly.parameters.window.fb_prestim = 20;
optoanaly.parameters.window.fb_stim = 20;
optoanaly.parameters.window.fb_step = 1;
optoanaly.parameters.window.fb_smooth = 1;

optoanaly.parameters.window.prestim = 20;
optoanaly.parameters.window.stim = 20;
optoanaly.parameters.window.step = 1;
optoanaly.parameters.window.smooth = 1;

%% Setting distance between tick windows
optoanaly.parameters.tick_dist = 10;

%% Make spiking timestamps centered at first stimulation of bout
for i = 1:optoanaly.timestamps.neuron_sum.num
  name = ['neuron_' num2str(i)];

  for j = 1:length(optoanaly.timestamps.opto_events.first)
    optoanaly.centered.fb.(name){j} = optoanaly.timestamps.(name).timestamps...
      (optoanaly.timestamps.(name).timestamps > ...
        (optoanaly.timestamps.opto_events.first(j)-optoanaly.parameters.window.fb_prestim) ...
      & optoanaly.timestamps.(name).timestamps < ...
        (optoanaly.timestamps.opto_events.first(j) + optoanaly.parameters.window.fb_stim))...
        - optoanaly.timestamps.opto_events.first(j);
  end
end

%% Make smooth psth of timestamps centered at first stimulation of bout
for i = 1:optoanaly.timestamps.neuron_sum.num
  name = ['neuron_' num2str(i)];
  optoanaly.parameters.window.fb_centers = -optoanaly.parameters.window.fb_prestim:optoanaly.parameters.window.fb_step:...
    optoanaly.parameters.window.fb_stim;

  for jj = 1:length(optoanaly.parameters.window.fb_centers)
    optoanaly.psth.fb.(name)(jj) = 0;
    for ii = 1:length(optoanaly.timestamps.opto_events.first)
      optoanaly.psth.fb.(name)(jj) = optoanaly.psth.fb.(name)(jj) + ((length(find(optoanaly.centered.fb.(name){ii} >...
        (optoanaly.parameters.window.fb_centers(jj)-(optoanaly.parameters.window.fb_smooth/2))...
        & optoanaly.centered.fb.(name){ii} < (optoanaly.parameters.window.fb_centers(jj)+...
        (optoanaly.parameters.window.fb_smooth/2)))))./((optoanaly.parameters.window.fb_smooth/1000)*...
        length(optoanaly.timestamps.opto_events.first)));
    end
  end
end


%% Make rasters centered at first stimulation of bout
for i = 1:optoanaly.timestamps.neuron_sum.num
  name = ['neuron_' num2str(i)];

  figure(i)
  subplot(3,2,1), hold on
  %hold on
  title('Rasters (first)')
  ylim([0 length(optoanaly.timestamps.opto_events.first)])
  ylabel('Trial')
  xlim([-optoanaly.parameters.window.fb_prestim ...
    optoanaly.parameters.window.fb_stim]);
  set(gca, 'Xtick', -optoanaly.parameters.window.fb_prestim:optoanaly.parameters.window.fb_stim)
 
  %making a cell for labeling ticks
  x_ind =  -optoanaly.parameters.window.fb_prestim:optoanaly.parameters.tick_dist:optoanaly.parameters.window.fb_stim;
  x_labels = cell(1,length(-optoanaly.parameters.window.fb_prestim:optoanaly.parameters.window.fb_stim));
  x = 1;
  for m = 1:optoanaly.parameters.tick_dist:length(x_labels)
    x_labels{m} = x_ind(x);
    x = x+1;
  end
  set(gca, 'XtickLabel', x_labels)

  xlabel('Time (msec)')

  %draw lines for when stimulus turned on
  line([0 0], get(gca, 'YLim'), 'Color', [0 0 1], 'Linewidth', 1, 'Linestyle', '-')

  %draw line for each spike in raster
  for k = 1:length(optoanaly.timestamps.opto_events.first)
    for j = 1:length(optoanaly.centered.fb.(name){k})
      line([optoanaly.centered.fb.(name){k}(j) optoanaly.centered.fb.(name){k}(j)], ...
        [k-1 k], 'LineWidth', 1, 'Color', [0 0 0])
    end
  end
  hold off
end

%% Make smooth PSTH figure centered at first stimulation of bout
for i = 1:optoanaly.timestamps.neuron_sum.num
    name = ['neuron_' num2str(i)];
    
    figure(i)
    subplot(3,2,3), hold on
    title('PSTH (first)')
    xlim([-optoanaly.parameters.window.fb_prestim ...
        optoanaly.parameters.window.fb_stim])
    set(gca, 'Xtick', -optoanaly.parameters.window.fb_prestim:optoanaly.parameters.window.fb_stim)
    
    %making a cell for labeling ticks
    x_ind =  -optoanaly.parameters.window.fb_prestim:optoanaly.parameters.tick_dist:optoanaly.parameters.window.fb_stim;
    x_labels = cell(1,length(-optoanaly.parameters.window.fb_prestim:optoanaly.parameters.window.fb_stim));
    x = 1;
    for m = 1:optoanaly.parameters.tick_dist:length(x_labels)
      x_labels{m} = x_ind(x);
      x = x+1;
    end
    set(gca, 'XtickLabel', x_labels)

    ylabel('Spikes/sec')
    xlabel('Time (msec)')
    
    %plot smooth psth
    plot(optoanaly.parameters.window.fb_centers, optoanaly.psth.fb.(name), 'LineWidth', 1, 'Color', [1 0 0])
    
    %draw line for when stimulus is turned on
    line([0 0], get(gca, 'YLim'), 'Color', [0 0 1], 'LineWidth',1, 'Linestyle', '--')
    hold off
end

%% Make spiking timestamps centered at every stimulation
for i = 1:optoanaly.timestamps.neuron_sum.num
  name = ['neuron_' num2str(i)];

  for j = 1:length(optoanaly.timestamps.opto_events.all)
    optoanaly.centered.all.(name){j} = optoanaly.timestamps.(name).timestamps...
      (optoanaly.timestamps.(name).timestamps > ...
        (optoanaly.timestamps.opto_events.all(j)-optoanaly.parameters.window.prestim) ...
      & optoanaly.timestamps.(name).timestamps < ...
        (optoanaly.timestamps.opto_events.all(j) + optoanaly.parameters.window.stim))...
        - optoanaly.timestamps.opto_events.all(j);
  end
end

%% Make smooth psth of timestamps centered at every stimulation
for i = 1:optoanaly.timestamps.neuron_sum.num
  name = ['neuron_' num2str(i)];
  optoanaly.parameters.window.centers = -optoanaly.parameters.window.prestim:optoanaly.parameters.window.step:...
    optoanaly.parameters.window.stim;

  for jj = 1:length(optoanaly.parameters.window.centers)
    optoanaly.psth.all.(name)(jj) = 0;
    for ii = 1:length(optoanaly.timestamps.opto_events.all)
      optoanaly.psth.all.(name)(jj) = optoanaly.psth.all.(name)(jj) + ((length(find(optoanaly.centered.all.(name){ii} >...
        (optoanaly.parameters.window.centers(jj)-(optoanaly.parameters.window.smooth/2))...
        & optoanaly.centered.all.(name){ii} < (optoanaly.parameters.window.centers(jj)+...
        (optoanaly.parameters.window.smooth/2)))))./((optoanaly.parameters.window.smooth/1000)*...
        length(optoanaly.timestamps.opto_events.all)));
    end
  end
end

%% Make rasters centered at every stimulation
for i = 1:optoanaly.timestamps.neuron_sum.num
  name = ['neuron_' num2str(i)];

  figure(i)
  subplot(3,2,2), hold on
  %hold on
  title(['Rasters (all)'])
  ylim([0 length(optoanaly.timestamps.opto_events.all)])
  ylabel('Trial')
  xlim([-optoanaly.parameters.window.prestim optoanaly.parameters.window.stim]);
  set(gca, 'Xtick', -optoanaly.parameters.window.prestim:optoanaly.parameters.window.stim)
  
  %making a cell for labeling ticks
  x_ind =  -optoanaly.parameters.window.prestim:optoanaly.parameters.tick_dist:optoanaly.parameters.window.stim;
  x_labels = cell(1,length(-optoanaly.parameters.window.prestim:optoanaly.parameters.window.stim));
  x = 1;
  for m = 1:optoanaly.parameters.tick_dist:length(x_labels)
    x_labels{m} = x_ind(x);
    x = x+1;
  end
  set(gca, 'XtickLabel', x_labels)
  xlabel('Time (msec)')

  %draw lines for when stimulus turned on
  line([0 0], get(gca, 'YLim'), 'Color', [0 0 1], 'Linewidth', 1, 'Linestyle', '-')

  %draw line for each spike in raster
  for k = 1:length(optoanaly.timestamps.opto_events.all)
    for j = 1:length(optoanaly.centered.all.(name){k})
      line([optoanaly.centered.all.(name){k}(j) optoanaly.centered.all.(name){k}(j)], ...
        [k-1 k], 'LineWidth', 1, 'Color', [0 0 0])
    end
  end
  hold off
end

%% Make smooth PSTH figure centered at every stimulation
for i = 1:optoanaly.timestamps.neuron_sum.num
    name = ['neuron_' num2str(i)];
    
    figure(i)
    subplot(3,2,4), hold on
    title('PSTH (first)')
    xlim([-optoanaly.parameters.window.prestim ...
        optoanaly.parameters.window.stim])
    set(gca, 'Xtick', -optoanaly.parameters.window.prestim:optoanaly.parameters.window.stim)

    %making a cell for labeling ticks
    x_ind =  -optoanaly.parameters.window.prestim:optoanaly.parameters.tick_dist:optoanaly.parameters.window.stim;
    x_labels = cell(1,length(-optoanaly.parameters.window.prestim:optoanaly.parameters.window.stim));
    x = 1;
    for m = 1:optoanaly.parameters.tick_dist:length(x_labels)
      x_labels{m} = x_ind(x);
      x = x+1;
    end
    set(gca, 'XtickLabel', x_labels)

    ylabel('Spikes/sec')
    xlabel('Time (msec)')
    
    %plot smooth psth
    plot(optoanaly.parameters.window.centers, optoanaly.psth.all.(name), 'LineWidth', 1, 'Color', [1 0 0])
    
    %draw line for when stimulus is turned on
    line([0 0], get(gca, 'YLim'), 'Color', [0 0 1], 'LineWidth',1, 'Linestyle', '--')
    hold off
end

%% Calculate ISIs and plot them
%Only plots ISIs less than the number specified
optoanaly.parameters.isi_maxplot = 120;

for i = 1:optoanaly.timestamps.neuron_sum.num
  name = ['neuron_' num2str(i)];

  optoanaly.isi.(name) = diff(optoanaly.timestamps.(name).timestamps);

  % Plot ISIs
  [N,X] = hist(optoanaly.isi.(name)...
    (optoanaly.isi.(name) < optoanaly.parameters.isi_maxplot), optoanaly.parameters.isi_maxplot);

  figure(i)
  subplot(3,2,5:6), hold on
  title('Interspike Intervals')
  ylabel('Proportion')
  xlabel('Time (msecs)')
  %set(gca, 'XTick', 0:optoanaly.parameters.isi_maxplot)
  %%making a cell for labeling ticks
  %x_ind =  0:optoanaly.parameters.tick_dist:optoanaly.parameters.isi_maxplot;  
  %x_labels = cell(1,length(0:optoanaly.parameters.isi_maxplot));
  %x = 1;
  %for m = 1:optoanaly.parameters.tick_dist:length(x_labels)
  %  x_labels{m} = x_ind(x);
  %  x = x+1;
  %end
  %set(gca, 'XtickLabel', x_labels)


  plot(X,N./sum(N))

  %save figure
  saveas(figure(gcf),['Optoanaly_of_' name],'fig')
  saveas(figure(gcf),['Optoanaly_of_' name],'epsc')
  saveas(figure(gcf),['Optoanaly_of_' name],'pdf')
end


