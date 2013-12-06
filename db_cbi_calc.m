function [cbi, ods_count, data, header] = db_cbi_calc(su_or_mu)
%db_cbi_calc Calculates CBI and distribution of ODS scores
%Must be in good_cells folder for an animal. Can specify
%if analysis should be done on single units, multi-units,
%or both. Saves ods_count and cbi to a .csv file and a figure
%comparing FR, cdf of ODI, and hist of ODS.
%
% example usage:  >>[cbi, ods_count, data, header] = db_cbi_calc('all')
%                 >>[cbi, ods_count, data, header]  = db_cbi_calc('su')
%                 >>[cbi, ods_count, data, header]  = db_cbi_calc('mu')
%
% CBI = [(n1-n7)+(2/3)(n2-n6)+(1/3)(n3-n5)+N]/2N
%
%WARNING! Uses 'annotation' function, which does not exist in octave
%so program will not run all the way through. Won't automatically
%save figure.
%
%Written by DM Brady Dec 2013

%% Checks to make sure you are in a 'good_cells' folder
curr_dir = pwd;

%if isempty(strfind(curr_dir,'good_cells')) == 1
%  display('You are not in a "good_cells" folder please change to the right directory');
%  display('Program terminated')
%  return
%end

%% Checks to make sure you have a good_cells.csv file to load
contents = dir([curr_dir '/good_cells.csv']);

if isempty(contents)
  display('You do not have a good_cells.csv file, have you run db_suanaly_v3 and saved the good cells?')
  display('Program terminated')
  return
end

%% Creates a word search based on whether you want to analyze SU, MU, or both
if strcmpi(su_or_mu, 'all')
  su_or_mu = 'u';
  appended = '_all';
elseif strcmpi(su_or_mu, 'su') || strcmpi(su_or_mu, 'mu')
  appended = ['_' su_or_mu];
else
  display('You did not enter "all", "su", or "mu". Please do so.')
  display('Program terminated')
end

%% Goes through the file and loads the relevant data
fid = fopen('good_cells.csv','r');
current_line = fgetl(fid);

% gives header information in a cell (does not include 'SU/MU')
commas = [0, strfind(current_line,',')];
for i = 1:length(commas)-1
  header{i} = current_line(commas(i)+1:commas(i+1)-1);
end

% goes through each line of .csv file and gets relevant info
current_line = fgetl(fid);
i = 1;
while ischar(current_line)
  if isempty(strfind(current_line,su_or_mu))
    current_line = fgetl(fid);
  else
    commas = [0, strfind(current_line,',')];
    for j = 1:length(commas)-1
      data(i,j) = str2num(current_line(commas(j)+1:commas(j+1)-1));
    end
    i = i+1;
    current_line = fgetl(fid);
  end
end
fclose(fid);

%% Plot normalized firing rates (contra vs. ipsi)
figure

subplot(2,2,1), hold on
title('Contra vs. Ipsi FR (Normalized)')
xlabel('Contra FR')
ylabel('Ipsi FR')
all_data = [data(:,strcmpi(header,'Contra_delta')); data(:,strcmpi(header,'Ipsi_delta'))];
if sum(all_data > 0) == length(all_data)
  xlim([0 1]);
  ylim([0 1]);
else
  xlim([-1 1]);
  ylim([-1 1]);
end
line(get(gca,'XLim'),get(gca,'YLim'),'Color',[1 0 0], 'LineWidth', 2,'LineStyle',':')
plot(data(:,strcmpi(header,'Contra_delta')), ...
  data(:,strcmpi(header,'Ipsi_delta')),'o','MarkerSize',7,'MarkerFaceColor',[0 0 1])
hold off


%% CDF plot of ODI scores
edges = -1:0.01:1;
count = histc(data(:,strcmpi(header,'ODI')),edges);
subplot(2,2,2), hold on
title('CDF of ODI')
xlabel('ODI score')
ylabel('Cumulative Fraction')
line([0 0], [0 1], 'Color', [1 0 0],'LineStyle',':', 'LineWidth',2)
plot(edges,cumsum(count)./sum(count),'LineWidth',2)
hold off

%% Histogram of ODS scores
ods_edges = 1:1:7;
ods_count = histc(data(:,strcmpi(header,'ODS')),ods_edges);
ods_count = ods_count';
subplot(2,2,3), hold on
title('ODS Histogram')
ylabel('Proportion')
xlabel('ODS Score')
bar(ods_edges, ods_count/sum(ods_count))
hold off


%% Calculates CBI score
% CBI = [(n1-n7)+(2/3)(n2-n6)+(1/3)(n3-n5)+N]/2N

cbi = [(ods_count(1)-ods_count(7))+...
       (2/3)*(ods_count(2)-ods_count(6))+...
       (1/3)*(ods_count(3)-ods_count(5))+...
       sum(ods_count)]/(2*sum(ods_count));

%% Writes ods_count and cbi to csv file
to_write = [ods_count, cbi];
csvwrite(['ODS_CBI' appended '.csv'],to_write);

annotation('textbox',[.6 .25 .3 .2],'FitBoxtoText','on','String',...
  {['    CBI: ' num2str(cbi)];...
  [];...
  ['    ODS  count: '];...
  [num2str(ods_edges) '  N'];...
  [num2str(ods_count) '   ' num2str(sum(ods_count))]});

saveas(figure(gcf),['Summary_CBI' appended],'fig')
saveas(figure(gcf),['Summary_CBI' appended],'pdf')
saveas(figure(gcf),['Summary_CBI' appended],'epsc')
