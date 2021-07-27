% TODO:
% squeeze/expand the timestamps of laser pulses to make them all a
% consistent length.
% Make is so that program determines automatically the laser on time
close all

addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\helper_functions'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\Spikes'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\npy-matlab'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\MatlabImportExport_v6.0.0'));

%% user inputs
% Neuralynx data path and session info:
%neuralynx_dir = 'F:\2020-05-16_08-07-57\12828547364_13295553145';
%timestamps_to_analyze = [12828547364 13031476251];
neuralynx_dir = 'F:\2020-05-12_11-56-00\31773437289_32247412068';
timestamps_to_analyze = [32104724156 32247412068];


PROBE_CONTACT_SPACING = 50; % set to 20 if H3 probe; 50 if L3 probe (um)
NLX_LFP_CHANNELS = 1:64;
CHANNEL_TO_READ_LFP = 1;
TIME_CONV_FACTOR = 10^-6; % 10^-6 means convert nlx units to seconds
VOLTAGE_CONV_FACTOR = 10^6; % 1 means output in volts, 1000 in mV, 10^6 in uV

LASER_PULSE_TIME = 0.1; % length of time laser was on for each pulse, s
%TODO:
%Calculate this from the data, it's not always the same...
PRE_LASER_TIME =  0.1; % length of time you want to look at before the start of laser, s
POST_LASER_TIME = 0.1; % length of time you want to look at after the stop of the laser, s
USE_MULTI_UNIT = true;
BIN_SIZE = 0.005; % bin size for binning mean firing rates, in seconds

numChannels = numel(NLX_LFP_CHANNELS);
kilosort_dir = fullfile(neuralynx_dir,'kilosort');
nlx_recording_times = strsplit(neuralynx_dir,'\'); nlx_recording_times = nlx_recording_times{end};
%% GET THE LFP TIMESTAMPS FROM THE NEURALYNX RECORDING
[lfp_ts, ~, sample_frequencies, number_of_valid_samples, lfp_samples, ~] = Nlx2MatCSC(fullfile(neuralynx_dir,strcat('CSC',num2str(CHANNEL_TO_READ_LFP),'_',nlx_recording_times,'.ncs')),[1 1 1 1 1],1,1,[]);
lfp_ts = lfp_ts.*TIME_CONV_FACTOR; %convert from us to s

[value, lfp_segment_start] = min(abs(lfp_ts - timestamps_to_analyze(1)*TIME_CONV_FACTOR));
[value, lfp_segment_end] = min(abs(lfp_ts - timestamps_to_analyze(2)*TIME_CONV_FACTOR));

lfp_ts = lfp_ts(lfp_segment_start:lfp_segment_end);
% redo the timestamps
lfp_ts_linearized = nan(512*size(lfp_ts,2),1);
for j = 1:numel(lfp_ts)
    if j < numel(lfp_ts)
        csc_dt = (lfp_ts(j+1) - lfp_ts(j))/512;
    else
        csc_dt = (lfp_ts(j) - lfp_ts(j-1))/512;
    end
    lfp_ts_linearized(512*(j-1)+1:512*j) = lfp_ts(j):csc_dt:lfp_ts(j)+csc_dt*511;
end


% % % get the sample frequence
% Fs = mode(sample_frequencies);
% % Fs = mean(1./(diff(lfp_ts)./512));
% 
% % extract information from header
% hdr = readCSCHeader(header);
% 
% % convert units
% lfp_samples = lfp_samples.*VOLTAGE_CONV_FACTOR.*hdr.ADBitVolts;
% 
% % construct within-block tvec
% nSamplesPerBlock = size(lfp_samples,1);
% block_tvec = (0:(1/Fs):((nSamplesPerBlock-1)/Fs));
% 
% % slow but tractable: loop over blocks remembering to initialize variables
% % first
% 
% lfp_ts_linearized = nan(numel(lfp_samples),1);
% 
% idx = 1; % move this along as we go through the loop over all samples
% 
% nBlocks = length(lfp_ts);
% bad_blocks = 0; % counter
% for iB = 1:nBlocks    
%     nvs = number_of_valid_samples(iB);
%     if nvs ~= 512, bad_blocks = bad_blocks + 1; end
%     current_data = lfp_samples(1:nvs,iB);
%     current_time = block_tvec(1:nvs)+lfp_ts(iB);    
%     lfp_ts_linearized(idx:idx+nvs-1) = current_time;   
%     idx = idx + nvs;   
% end % of block loop
session_length = (max(lfp_ts_linearized) - min(lfp_ts_linearized)); % s

%% DETERMINE WHEN THE LASER WAS ON
%% LOAD THE EVENTS FILE
neuralynx_event_file = fullfile(neuralynx_dir,strcat('Events_',nlx_recording_times,'.csv'));
events = readtable(neuralynx_event_file); 
event_segment_start = find(events.TimeStamp == timestamps_to_analyze(1));
event_segment_end = find(events.TimeStamp == timestamps_to_analyze(2));
events = events(event_segment_start:event_segment_end,:);
ts_laser_on_start = events.TimeStamp(events.EventID == 11 & events.TTLValue == 1).*TIME_CONV_FACTOR;
ts_laser_off_start = events.TimeStamp(events.EventID == 11 & events.TTLValue == 0).*TIME_CONV_FACTOR;

% if the laser was on at the beginning, get rid of the first "off" pulse
if ts_laser_off_start(1) < ts_laser_on_start(1) 
    ts_laser_off_start(1) = [];
end

% if the laser ended by turning on, get rid of the last laser on pulse
if length(ts_laser_on_start) > length(ts_laser_off_start)
    ts_laser_on_start(end) = [];
end

numLaserPulses = numel(ts_laser_on_start);
%% Find the indices associated with each segmant (pre-laser pulse, laser pulse, post-laser pulse)
ind_laser_on_start = nan(numLaserPulses,1);
ind_post_laser_start = nan(numLaserPulses,1);
ind_pre_laser_start = nan(numLaserPulses,1);
ind_post_laser_stop = nan(numLaserPulses,1);

for pulse = 1:numLaserPulses
    % first idx of laser on
    [~,ind_laser_on_start(pulse)] = min(abs(lfp_ts_linearized - ts_laser_on_start(pulse)));
    % first idx following offset of laser
    [~,ind_post_laser_start(pulse)] = min(abs(lfp_ts_linearized - ts_laser_off_start(pulse)));
    % idx for the beginning of the laser_off period prior to a laser pulse
    % OLD Way:
    time_pre_laser_start = lfp_ts_linearized(ind_laser_on_start(pulse)) - PRE_LASER_TIME;
    [~, ind_pre_laser_start(pulse)] = min(abs(lfp_ts_linearized - time_pre_laser_start));
    % idx for the end of the laser_off period following a laser pulse
    time_post_laser_stop = lfp_ts_linearized(ind_post_laser_start(pulse)) + POST_LASER_TIME;
    [~, ind_post_laser_stop(pulse)] = min(abs(lfp_ts_linearized - time_post_laser_stop));    
end

ind_laser_on_stop = ind_post_laser_start - 1;
ind_pre_laser_stop = ind_laser_on_start - 1;

if ind_pre_laser_start(1) < 1
    ind_laser_on_start(1) = [];
    ind_post_laser_start(1) = [];
    ind_pre_laser_start(1) = [];
    ind_post_laser_stop(1) = [];
    ind_laser_on_stop(1) = [];
    ind_pre_laser_stop(1) = [];
end

if ind_post_laser_stop(end) > length(lfp_ts_linearized)
    ind_laser_on_start(end) = [];
    ind_post_laser_start(end) = [];
    ind_pre_laser_start(end) = [];
    ind_post_laser_stop(end) = [];
    ind_laser_on_stop(end) = [];
    ind_pre_laser_stop(end) = [];
end
    

%% Find the timestamps associated with each segmant (pre-laser pulse, laser pulse, post-laser pulse)
ts_laser_on_start = lfp_ts_linearized(ind_laser_on_start);
ts_post_laser_start = lfp_ts_linearized(ind_post_laser_start);
ts_pre_laser_start = lfp_ts_linearized(ind_pre_laser_start);
ts_post_laser_stop = lfp_ts_linearized(ind_post_laser_stop);
ts_laser_on_stop = lfp_ts_linearized(ind_laser_on_stop);
ts_pre_laser_stop = lfp_ts_linearized(ind_pre_laser_stop);

numLaserPulses = numel(ts_laser_on_start);

%%
% load the silicon probe cluster_info file containing basic info about each
% cluster
% remove noise clusters
[~,header,raw] = tsvread(fullfile(kilosort_dir,'cluster_info.tsv'));
noise_clusters = find(strcmp(raw(2:end,6),'noise'));
depth = raw(2:end,4);
depth = cellfun(@str2num, depth);
depth(noise_clusters) = [];

% load the phy output containing timestamps for each cluster
sp = loadKSdir(kilosort_dir);
spike_times = sp.st; %
spike_ids = sp.clu;
spike_inds = round(spike_times*32000);


% spike_inds = spike_inds(spike_inds >= (lfp_segment_start - 1)*512 + 1 & spike_inds <= lfp_segment_end*512);
% spike_ids = spike_ids(spike_inds >= (lfp_segment_start - 1)*512 + 1 & spike_inds <= lfp_segment_end*512);
spike_inds = spike_inds(spike_inds >= (lfp_segment_start - 1)*512 & spike_inds <= lfp_segment_end*512);
spike_ids = spike_ids(spike_inds >= (lfp_segment_start - 1)*512 & spike_inds <= lfp_segment_end*512);

% Trying to figure out timing issue
% We only want the spikes from the segment that we are currently analyzing:
% 
%spike_inds = spike_inds - ((lfp_segment_start - 1)*512 + 1);
spike_inds = spike_inds - ((lfp_segment_start - 1)*512 );
spike_times = lfp_ts_linearized(spike_inds);


cluster_id = sp.cids;
cluster_type = sp.cgs;

% pull out the cells or multi-units

if USE_MULTI_UNIT == 1
    good_cluster_id = cluster_id(cluster_type==1 | cluster_type==2);
    good_cluster_depth = depth(cluster_type==1 | cluster_type==2);
else
    good_cluster_id = cluster_id(cluster_type==2);
    good_cluster_depth = depth(cluster_type==2);
end

% sort by depth
[~, sorted_depth_idx] = sort(good_cluster_depth);
numCells = numel(good_cluster_id);
%% Make a raster plot of spikes at the time before, during, and after each laser pulse

numBins = round((PRE_LASER_TIME + LASER_PULSE_TIME + POST_LASER_TIME)/BIN_SIZE);
binEdges = -PRE_LASER_TIME:BIN_SIZE:LASER_PULSE_TIME + POST_LASER_TIME;
binCenters = (binEdges(1:end-1)+binEdges(2:end))/2;
smoothing_window_length = 1; 
w = gausswin(smoothing_window_length./BIN_SIZE)/sum(gausswin(smoothing_window_length/BIN_SIZE));


pre_laser_spikes = cell(numCells,numLaserPulses);
laser_spikes = cell(numCells,numLaserPulses);
post_laser_spikes = cell(numCells,numLaserPulses);
rel_pre_laser_spikes = cell(numCells,numLaserPulses);
rel_laser_spikes = cell(numCells,numLaserPulses);
rel_post_laser_spikes = cell(numCells,numLaserPulses);
all_near_laser_spikes = cell(numCells, numLaserPulses);
all_near_laser_spikes_hist = cell(numCells,1);
num_spikes_no_laser = nan(numCells,numLaserPulses);
num_spikes_laser = nan(numCells,numLaserPulses);

for this_cell = 1:numCells
    cell_spike_times = spike_times(spike_ids == good_cluster_id(this_cell));
    cell_spike_hist = nan(numLaserPulses,numBins);
    
    for pulse = 1:numLaserPulses
        pre_laser_start = ts_pre_laser_start(pulse);
        pre_laser_stop = ts_pre_laser_stop(pulse);
        laser_start = ts_laser_on_start(pulse);
        laser_stop = ts_laser_on_stop(pulse);
        post_laser_start = ts_post_laser_start(pulse);
        post_laser_stop = ts_post_laser_stop(pulse);
        
        pre_laser_spikes{this_cell,pulse} = cell_spike_times(cell_spike_times >= pre_laser_start & cell_spike_times <= pre_laser_stop);
        laser_spikes{this_cell,pulse} = cell_spike_times(cell_spike_times >= laser_start & cell_spike_times <= laser_stop);
        post_laser_spikes{this_cell,pulse} = cell_spike_times(cell_spike_times >= post_laser_start & cell_spike_times <= post_laser_stop);
        
        num_spikes_no_laser(this_cell,pulse) = numel(pre_laser_spikes{this_cell,pulse}) + numel(post_laser_spikes{this_cell,pulse});
        num_spikes_laser(this_cell,pulse) = numel(laser_spikes{this_cell,pulse});
        
        rel_pre_laser_spikes = (cell_spike_times(cell_spike_times >= pre_laser_start & cell_spike_times <= pre_laser_stop) - laser_start);
        rel_laser_spikes = (cell_spike_times(cell_spike_times >= laser_start & cell_spike_times <= laser_stop) - laser_start);
        rel_post_laser_spikes = (cell_spike_times(cell_spike_times >= post_laser_start & cell_spike_times <= post_laser_stop) - laser_start);
        
        all_near_laser_spikes{this_cell,pulse} = [rel_pre_laser_spikes' rel_laser_spikes' rel_post_laser_spikes'];
        spike_hist = hist(all_near_laser_spikes{this_cell,pulse},binCenters);      
        cell_spike_hist(pulse,:) = spike_hist;
    end
    
    all_near_laser_spikes_hist{this_cell} = cell_spike_hist;
end


mean_spike_hist = cell2mat(cellfun(@mean, all_near_laser_spikes_hist, 'UniformOutput', false));
std_spike_hist = cell2mat(cellfun(@std, all_near_laser_spikes_hist, 'UniformOutput',false));
% bar(mean_spike_hist{1}) will plot the first cell

mean_firing_rate_hist = mean_spike_hist./BIN_SIZE;
std_firing_rate_hist = std_spike_hist./BIN_SIZE;

% smoothing probably won't work here, because it causes a bit of a time
% shift... but if you do want to smooth with a Gaussian window do the
% following:
smoothed_firing_rates = cell(numCells,1);
for this_cell = 1:numCells
    this_cell_smoothed_firing_rates = nan(numLaserPulses,numBins);
    for pulse = 1:numLaserPulses
        this_cell_smoothed_firing_rates(pulse,:) = filter(w,1,all_near_laser_spikes_hist{this_cell}(pulse,:));
    end
    smoothed_firing_rates{this_cell} = this_cell_smoothed_firing_rates;
end

% This does't work: I was trying to normalize by the pre-laser FR, but if that's very low (like 0.01), 
% it returns a huge value.
mean_firing_rate_laser_off = mean(mean_spike_hist(:,binCenters < 0),2)./BIN_SIZE;
mean_firing_rate_laser_on = mean(mean_spike_hist(:,binCenters > 0 & binCenters < LASER_PULSE_TIME),2)./BIN_SIZE;
% mean_firing_rate_hist_normalized = mean_spike_hist./(mean(mean_spike_hist(:,binCenters < 0),2));

%mean_firing_rate_hist_normalized = mean_spike_hist./(mean(mean_firing_rate_laser_off + mean_firing_rate_laser_on));
mean_firing_rate_hist_normalized = mean_firing_rate_hist./(mean_firing_rate_laser_off);


% for this_cell = 1:numCells
% figure()
% subplot(1,3,1)
% bar(binCenters,spike_histograms(this_cell,:));
% xlabel('Time (s) from Laser onset')
% ylabel('Spikes/s')
% 
% subplot(1,3,2)
% bar(binCenters,single_lap_spike_histograms(this_cell,:));
% xlabel('Time (s) from Laser onset')
% ylabel('# spikes')
% 
% y = smoothed_firing_rates(this_cell,:); % your mean vector;
% x = binCenters;
% std_dev = std_spike_hist;
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% subplot(1,3,3)
% fill(x2, inBetween, [0.6,0.6,0.6]);
% hold on;
% plot(x, y, 'k', 'LineWidth', 2);
% xlabel('Time (s) from Laser onset')
% ylabel('Smoothed firing rate')
% end

sorted_depth = depth(sorted_depth_idx);
all_near_laser_spikes_sorted_by_depth = all_near_laser_spikes(sorted_depth_idx,:);
all_near_laser_spikes_hist_sorted_by_depth = all_near_laser_spikes_hist(sorted_depth_idx,:);
firing_rate_laser_off = num_spikes_no_laser./(PRE_LASER_TIME+POST_LASER_TIME);
firing_rate_laser_on = num_spikes_laser./(LASER_PULSE_TIME);
mean_firing_rate_laser_off = mean(firing_rate_laser_off,2);
mean_firing_rate_laser_on = mean(num_spikes_laser,2)./(LASER_PULSE_TIME);
modulation_index = (mean_firing_rate_laser_off - mean_firing_rate_laser_on)./(mean_firing_rate_laser_off + mean_firing_rate_laser_on);
modulation_index_sorted_by_depth = modulation_index(sorted_depth_idx);

% Make a histogram of modulaton index for all cells
figure()
bin_centers = linspace(-0.9,0.9,10);
hist(modulation_index,bin_centers)
xlim([-1 1])
title('Modulation Index')
ylabel('Count')
xlabel('(FR Laser ON - FR Laser OFF)/(FR Laser ON + FR Laser OFF)')
saveas(gcf,fullfile(kilosort_dir,'Modulation Index'),'png')

modulation_bins = (-1:0.1:1);
modulation_value = nan(numCells,1);
for this_cell = 1:numCells
    for bin = 1:length(modulation_bins) - 1
        if bin < length(modulation_bins) - 1
            if modulation_index(this_cell) >= modulation_bins(bin) && modulation_index(this_cell) < modulation_bins(bin + 1)
                modulation_value(this_cell) = bin;
            end
        elseif bin == length(modulation_bins) - 1
            if modulation_index(this_cell) >= modulation_bins(bin) && modulation_index(this_cell) <= modulation_bins(bin + 1)
                modulation_value(this_cell) = bin;
            end
        end
        
    end
end
modulation_value_sorted_by_depth = modulation_value(sorted_depth_idx);

% Determine which cells are significantly impacted by the Laser
significant_cells = zeros(numCells,1);
for this_cell = 1:numCells
    pval = signrank(firing_rate_laser_off(this_cell,:),firing_rate_laser_on(this_cell,:));
    if pval < 0.05
        significant_cells(this_cell) = 1;
    end
end
%%
figure()
significantly_inhibited_cells = find(modulation_index > 0 & significant_cells == 1);
mean_response_inhibited_cells = mean(mean_firing_rate_hist_normalized(significantly_inhibited_cells,:));
std_dev = std(mean_firing_rate_hist_normalized(significantly_inhibited_cells,:));
std_err = std_dev/sqrt(length(significantly_inhibited_cells));

y = mean_response_inhibited_cells; % your mean vector;
x = binCenters;
fill([0, mean(ts_laser_on_stop - ts_laser_on_start),  mean(ts_laser_on_stop - ts_laser_on_start), 0 ],[1.9 1.9 2 2],'r', 'LineStyle','none')
hold on
curve1 = y + std_err;
curve2 = y - std_err;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.6,0.6,0.6]);
hold on;
plot(x, y, 'k', 'LineWidth', 2);
hold on
xlabel('Time from laser onset (s)','FontSize',15)
ylabel('Mean response of significantly inhibited cells, normalized to pre-laser firing rate','FontSize',15)
xticks(-1*PRE_LASER_TIME:LASER_PULSE_TIME/5:LASER_PULSE_TIME+POST_LASER_TIME)
xlim([-1*PRE_LASER_TIME LASER_PULSE_TIME+POST_LASER_TIME]) 
saveas(gcf,fullfile(kilosort_dir,['Mean_response_inhibited_cells_',num2str(timestamps_to_analyze(1)),'_',num2str(timestamps_to_analyze(2))]),'png')

%%
figure()
significantly_excited_cells = find(modulation_index < 0 & significant_cells == 1);
mean_response_excited_cells = mean(mean_firing_rate_hist_normalized(significantly_excited_cells,:),1);
std_dev = std(mean_firing_rate_hist_normalized(significantly_excited_cells,:),[],1);
std_err = std_dev/sqrt(length(significantly_excited_cells));

y = mean_response_excited_cells; % your mean vector;
x = binCenters;
fill([0, mean(ts_laser_on_stop - ts_laser_on_start),  mean(ts_laser_on_stop - ts_laser_on_start), 0 ],[0 0 2 2],'r', 'LineStyle','none')
hold on
curve1 = y + std_err;
curve2 = y - std_err;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.6,0.6,0.6]);
hold on;
plot(x, y, 'k', 'LineWidth', 2);
xlabel('time (s) from Laser onset')
ylabel('Response, normalized')
xticks(-1*PRE_LASER_TIME:LASER_PULSE_TIME/10:LASER_PULSE_TIME+POST_LASER_TIME)
xlim([-1*PRE_LASER_TIME LASER_PULSE_TIME+POST_LASER_TIME]) 
saveas(gcf,fullfile(kilosort_dir,['Mean_response_excited_cells_',num2str(timestamps_to_analyze(1)),'_',num2str(timestamps_to_analyze(2))]),'png')
%%
significant_cells_sorted_by_depth = significant_cells(sorted_depth_idx);
percent_significant_cells = sum(significant_cells_sorted_by_depth)/length(significant_cells_sorted_by_depth)*100;

% Make a figure plotting the modulation index as a function of depth along
% the probe
figure()
significance_fill_color = nan(numCells,3);
for i = 1:numCells
    if(significant_cells_sorted_by_depth(i)) == 1
        significance_fill_color(i,:) = [1 0 0];
    else
        significance_fill_color(i,:) = [.6 .6 .6];
    end
end
s = scatter((sorted_depth-1).*PROBE_CONTACT_SPACING - 3150, modulation_index_sorted_by_depth, [], significance_fill_color);
s.MarkerFaceColor = [.6 .6 .6];
s.LineWidth = 1;
P = polyfit((sorted_depth-1).*PROBE_CONTACT_SPACING - 3150 ,modulation_index(sorted_depth_idx),1);
x0 = min((sorted_depth-1).*PROBE_CONTACT_SPACING - 3150) ; x1 = max((sorted_depth-1).*PROBE_CONTACT_SPACING - 3150) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
plot(xi,yi,'k','LineWidth',3);
ylabel('modulation index; +1 = inhibition, -1 = excitation','FontSize',15)
xlabel('distance from probe tip, uM','FontSize',15)
ylim([-1 1])
formatSpec = '%2.1f%% of cells were significantly modulated by light' ;
A1 = percent_significant_cells;
title_str = sprintf(formatSpec,A1);
title(title_str,'FontSize',15)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [.1, 0, .6, 1]);
saveas(gcf,fullfile(kilosort_dir,['Modulation Index versus Depth_',num2str(timestamps_to_analyze(1)),'_',num2str(timestamps_to_analyze(2))]),'png')


% Make a scatterplot showing spikes before, during, and after Laser Period
cmap = flipud(jet(length(modulation_bins)-1));
figure()

for this_cell = 1:numCells
    % pull out the cell containing the spikes on each trial for a given cell 
    spikes = all_near_laser_spikes_sorted_by_depth(this_cell,:);
    [xPoints, yPoints] = plotSpikeRaster(spikes,'PlotType','scatter');
    y_shift = size(all_near_laser_spikes_sorted_by_depth,2)*(numCells-this_cell);
    yPoints = yPoints + y_shift;
    this_cell_color = cmap(modulation_value_sorted_by_depth(this_cell),:);
    plot(xPoints,yPoints,'.k','MarkerSize',12,'color',this_cell_color)
    hold on
end

ylim([0 numLaserPulses*numCells + 1])
title('Inhibition of MEC by Jaws','FontSize',15)
xticks(-1*PRE_LASER_TIME:LASER_PULSE_TIME/5:LASER_PULSE_TIME+POST_LASER_TIME)
xlabel('Time, relative to laser onset (s)','FontSize',15)
yticks(numLaserPulses/2:numLaserPulses:numLaserPulses*numCells-numLaserPulses/2)
yticklabels(numCells:-1:1)
ylabel('Cell number','FontSize',15)
axis tight
colormap(flipud(jet))
cb = colorbar;
%set(cb,'YDir','reverse')

ylabel(cb,'Modulation Index (FR laser ON - FR laser OFF)/(FR laser ON + FR laser OFF)','FontSize',15)
caxis([-1 1])
axl = gca;

yyaxis right
axr = gca;
axr.YLim = [0 numLaserPulses*numCells];
axr.YTick = numLaserPulses/2:numLaserPulses:numLaserPulses*numCells-numLaserPulses/2;
distance_along_probe = 3.125 - ((sorted_depth)./64).*3.125;
axr.YTickLabels = round(flipud(distance_along_probe),1);

ylabel('Distance from probe tip (mm)','FontSize',15)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [.2, 0, .5, 1]);
linkaxes([axl,axr],'y');
xlim([-1*PRE_LASER_TIME LASER_PULSE_TIME+POST_LASER_TIME]) 
set(gcf, 'renderer', 'painters')
keyboard
saveas(gcf,fullfile(kilosort_dir,['Raster Plot_',num2str(timestamps_to_analyze(1)),'_',num2str(timestamps_to_analyze(2))]),'png')

keyboard

% for this_cell = 1:numCells
%     for pulse = 1:numLaserPulses
% %         y_value = numLaserPulses*(numChannels - sorted_depth(this_cell)) + pulse;
%         y_value = numLaserPulses*(numCells + 1 - this_cell) + pulse - numLaserPulses;
%         scatter(bins,all_near_laser_spikes_hist_sorted_by_depth{this_cell,pulse}*y_value,4,cmap(sorted_modulation_value(this_cell),:),'filled')
%         hold on
%     end
%    
%    
% end

% xlim([min(bins) max(bins)]) 
% ylim([0 numLaserPulses*numCells])
% title('Inhibition of MEC by Jaws')
% xticks(-1*PRE_LASER_TIME:binSize*10:3*LASER_PULSE_TIME)
% xlabel('time, relative to laser onset (s)')
% yticks(numLaserPulses/2:numLaserPulses:numLaserPulses*numCells-numLaserPulses/2)
% yticklabels(numCells:-1:1)
% ylabel('cell number')
% axis tight

% colormap(jet)
% cb = colorbar;
% set(cb,'YDir','reverse')
% ylabel(cb,'Modulation Index')
% caxis([-1 1])
% 
% yyaxis right
% axr = gca;
% axr.YLim = [0 numLaserPulses*numCells];
% axr.YTick = numLaserPulses/2:numLaserPulses:numLaserPulses*numCells-numLaserPulses/2;
% axr.YTickLabels = flipud(sorted_depth);
% ylabel('Position Along Probe')
% saveas(gcf,fullfile(kilosort_dir,'Raster Plot'),'png')


