% analyzeProbeData
% cids = cluster group; 1 = MUA, 2 = Good, 3 = Unsorted
% sp = spike times
% sp.clu = clusterID
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\helper_functions'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\Spikes'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\npy-matlab'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\MatlabImportExport_v6.0.0'));

%% user inputs
% Neuralynx data path and session info:
neuralynx_dir = 'F:\2020-02-24_13-25-55\19025800573_19442134375';
NLX_LFP_CHANNELS = 1:64;
CHANNEL_TO_READ_LFP = 1;
TIME_CONV_FACTOR = 1; % 10^-6 means convert nlx units to seconds
VOLTAGE_CONV_FACTOR = 10^6; % 1 means output in volts, 1000 in mV, 10^6 in uV
PRE_LASER_TIME = 0.1*10^6; % length of time you want to look at before the start of laser, us
LASER_PULSE_TIME = 0.1*10^6; % length of time laser was on for each pulse, us
POST_LASER_TIME = 0.1*10^6; % length of time you want to look at after the stop of the laser, us
USE_MULTI_UNIT = true;


numChannels = numel(NLX_LFP_CHANNELS);
kilosort_dir = fullfile(neuralynx_dir,'kilosort');
nlx_recording_times = strsplit(neuralynx_dir,'\'); nlx_recording_times = nlx_recording_times{end};
%% GET THE LFP TIMESTAMPS FROM THE NEURALYNX RECORDING

[lfp_ts, ~, sample_frequencies, number_of_valid_samples, lfp_samples, header] = Nlx2MatCSC(fullfile(neuralynx_dir,strcat('CSC',num2str(CHANNEL_TO_READ_LFP),'_',nlx_recording_times,'.ncs')),[1 1 1 1 1],1,1,[]);

% get the sample frequence
Fs = mode(sample_frequencies);

% extract information from header
hdr = readCSCHeader(header);

% convert units
lfp_samples = lfp_samples .*VOLTAGE_CONV_FACTOR .*hdr.ADBitVolts;

% construct within-block tvec
nSamplesPerBlock = size(lfp_samples,1);
block_tvec = (0:(1/Fs):((nSamplesPerBlock-1)/Fs)).*10^6;

% slow but tractable: loop over blocks remembering to initialize variables
% first

lfp_ts_linearized = nan(numel(lfp_samples),1);

idx = 1; % move this along as we go through the loop over all samples

nBlocks = length(lfp_ts);
bad_blocks = 0; % counter
for iB = 1:nBlocks
    
    nvs = number_of_valid_samples(iB);
    if nvs ~= 512, bad_blocks = bad_blocks + 1; end
    
    current_data = lfp_samples(1:nvs,iB);
    current_time = block_tvec(1:nvs)+lfp_ts(iB);
    
    lfp_ts_linearized(idx:idx+nvs-1) = current_time;
    
    idx = idx + nvs;
    
end % of block loop
session_length = (max(lfp_ts_linearized) - min(lfp_ts_linearized))/(10^6); % s
%% LOAD ALL LFP DATA INTO A BIG MATRIX -- at some point make this a function

% all_lfp_samples_linearized = nan(numChannels,length(lfp_ts_linearized));
% % 
% for c = 1:numChannels
%     % load a neuralynx eeg file to get the timestamps
%     [lfp_ts, ~, sample_frequencies, number_of_valid_samples, lfp_samples, header] = Nlx2MatCSC(fullfile(neuralynx_dir,strcat('CSC',num2str(NLX_LFP_CHANNELS(c)),'_',nlx_recording_times,'.ncs')),[1 1 1 1 1],1,1,[]);
% 
%     % extract information from header
%     hdr = readCSCHeader(header);
%     
%     % convert units
%     lfp_samples = lfp_samples .*VOLTAGE_CONV_FACTOR.*hdr.ADBitVolts;
%     
%     % construct within-block tvec
%     nSamplesPerBlock = size(lfp_samples,1);
%     block_tvec = (0:(1/Fs):((nSamplesPerBlock-1)/Fs)).*10^6;
%     
%     % slow but tractable: loop over blocks remembering to initialize variables
%     % first
%     
%     lfp_samples_linearized = nan(numel(lfp_samples),1); % allocate memory and then fill; can trim later
%     
%     idx = 1; % move this along as we go through the loop over all samples
%     
%     nBlocks = length(lfp_ts);
%     bad_blocks = 0; % counter
%     for iB = 1:nBlocks
%         
%         nvs = number_of_valid_samples(iB);
%         if nvs ~= 512, bad_blocks = bad_blocks + 1; end
%         
%         current_data = lfp_samples(1:nvs,iB);
%         current_time = block_tvec(1:nvs)+lfp_ts(iB);
%         
%         lfp_samples_linearized(idx:idx+nvs-1) = current_data;
%         idx = idx + nvs;
%         
%     end % of block loop
%     all_lfp_samples_linearized(c,:) = lfp_samples_linearized;
% end

%% DETERMINE WHEN THE LASER WAS ON
%% LOAD THE EVENTS FILE
neuralynx_event_file = fullfile(neuralynx_dir,strcat('Events_',nlx_recording_times,'.csv'));
events = readtable(neuralynx_event_file); ts_laser_on_start = events.TimeStamp(events.EventID == 11 & events.TTLValue == 1);
ts_laser_off_start = events.TimeStamp(events.EventID == 11 & events.TTLValue == 0);
numLaserPulses = numel(ts_laser_on_start);

% if the laser ended by turning on, get rid of the last laser on pulse
if ts_laser_on_start(1) < ts_laser_off_start(1) && length(ts_laser_on_start) > length(ts_laser_off_start)
    ts_laser_on_start(end) = [];
    numLaserPulses = numLaserPulses - 1;
end
    
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
    ind_pre_laser_start(pulse) = ind_laser_on_start(pulse) - PRE_LASER_TIME*Fs/10^6;
    % idx for the end of the laser_off period following a laser pulse
    ind_post_laser_stop(pulse) = ind_post_laser_start(pulse) + POST_LASER_TIME*Fs/10^6 - 1;
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
%% Look at the LFP before, during, and after each laser pulse
% chunk_time = (PRE_LASER_TIME + LASER_PULSE_TIME + POST_LASER_TIME)/10^6; %s
% lfp_near_laser = nan(length(NLX_LFP_CHANNELS),Fs*chunk_time,numel(ts_laser_on_start));
% lfp_near_laser_normalized = nan(length(NLX_LFP_CHANNELS),Fs*chunk_time,numel(ts_laser_on_start));
% 
% for pulse = 1:numLaserPulses
%     
%     lfp_during_laser = all_lfp_samples_linearized(:,ind_laser_on_start(pulse):ind_laser_on_stop(pulse));
%     
%     % this is supposed to be 160000 samples (5 sec) but the actual laser ON times are inconsisent.
%     % Interpolate between the start and stop of the laser.
%     
%     x_samples = 1:size(lfp_during_laser,2);
%     xq_samples = linspace(1,size(lfp_during_laser,2),LASER_PULSE_TIME/(10^6)*Fs);
%     
%     lfp_during_laser_interp = nan(numChannels,length(xq_samples));
%     for channel = 1:size(lfp_during_laser,1)
%         lfp_during_laser_interp(channel,:) = interp1(x_samples,lfp_during_laser(channel,:),xq_samples);
%     end
%         
%     lfp_before_laser = all_lfp_samples_linearized(:,ind_pre_laser_start(pulse):ind_pre_laser_stop(pulse));
%     lfp_after_laser = all_lfp_samples_linearized(:,ind_post_laser_start(pulse):ind_post_laser_stop(pulse));
%     
%     lfp_near_laser(:,:,pulse) = [lfp_before_laser lfp_during_laser_interp lfp_after_laser];
%     lfp_near_laser_this_pulse = lfp_near_laser(:,:,pulse);
%     lfp_near_laser_normalized(:,:,pulse) = (lfp_near_laser_this_pulse-(min(lfp_near_laser_this_pulse,[],2)))./(max(lfp_near_laser_this_pulse,[],2)-min(lfp_near_laser_this_pulse,[],2));
% end
% mean_lfp_near_laser= mean(lfp_near_laser,3);
% mean_lfp_near_laser_normalized = mean(lfp_near_laser_normalized,3);
% 
% plot_channels = 1:9:65;
% figure
% for c = 1:length(plot_channels)
% subplot(length(plot_channels),1,c)
% hold on
% plot([PRE_LASER_TIME*Fs/10^6 PRE_LASER_TIME*Fs/10^6],[min(min(mean_lfp_near_laser_normalized(plot_channels,:), [], 2)) max(max(mean_lfp_near_laser_normalized(plot_channels,:), [], 2))],'-k')
% hold on
% plot([(PRE_LASER_TIME*Fs/10^6 + LASER_PULSE_TIME*Fs/10^6) (PRE_LASER_TIME*Fs/10^6 + LASER_PULSE_TIME*Fs/10^6)],[min(min(mean_lfp_near_laser_normalized(plot_channels,:), [], 2)) max(max(mean_lfp_near_laser_normalized(plot_channels,:), [], 2))],'-k')
% hold on
% plot(mean_lfp_near_laser_normalized(plot_channels(c),:),'r','LineWidth',2)
% ylabel(['channel ' num2str(plot_channels(c))])
% linkaxes
% end
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
spike_times = sp.st.*10^6; %change to microsec to be consistent with neuralynx
spike_times = spike_times + lfp_ts_linearized(1);

% pull out the cells or multi-units
cluster_id = sp.cids;

if USE_MULTI_UNIT == 1
    good_cluster_id = cluster_id(sp.cgs==1 | sp.cgs ==2);
    good_cluster_depth = depth(sp.cgs==1 | sp.cgs ==2);
else
    good_cluster_id = cluster_id(sp.cgs ==2);
    good_cluster_depth = depth(sp.cgs ==2);
end

% sort by depth
[~, sorted_depth_idx] = sort(good_cluster_depth);
numCells = numel(good_cluster_id);
%% Make a raster plot of spikes at the time before, during, and after each laser pulse
pre_laser_spikes = cell(numCells,numLaserPulses);
laser_spikes = cell(numCells,numLaserPulses);
post_laser_spikes = cell(numCells,numLaserPulses);
rel_pre_laser_spikes = cell(numCells,numLaserPulses);
rel_laser_spikes = cell(numCells,numLaserPulses);
rel_post_laser_spikes = cell(numCells,numLaserPulses);
all_near_laser_spikes = cell(numCells, numLaserPulses);
num_spikes_no_laser = nan(numCells,numLaserPulses);
num_spikes_laser = nan(numCells,numLaserPulses);
all_near_laser_spikes_hist = cell(numCells, numLaserPulses);
bins = -PRE_LASER_TIME/10^6:0.005:PRE_LASER_TIME/10^6 + LASER_PULSE_TIME/10^6 + POST_LASER_TIME/10^6;

for this_cell = 1:numCells
    cell_spike_times = spike_times(sp.clu == good_cluster_id(this_cell));
    for pulse = 1:numLaserPulses
        pre_laser_start = ts_pre_laser_start(pulse);
        pre_laser_stop = ts_pre_laser_stop(pulse);
        laser_start = ts_laser_on_start(pulse);
        laser_stop = ts_laser_on_stop(pulse);
        post_laser_start = ts_post_laser_start(pulse);
        post_laser_stop = ts_post_laser_stop(pulse);
        
        rel_pre_laser_spikes{this_cell,pulse} = (cell_spike_times(cell_spike_times >= pre_laser_start & cell_spike_times <= pre_laser_stop) - laser_start)/(10^6);
        rel_laser_spikes{this_cell,pulse} = (cell_spike_times(cell_spike_times >= laser_start & cell_spike_times <= laser_stop) - laser_start)/(10^6);
        rel_post_laser_spikes{this_cell,pulse} = (cell_spike_times(cell_spike_times >= post_laser_start & cell_spike_times <= post_laser_stop) - laser_start)/(10^6);
        
        pre_laser_spikes{this_cell,pulse} = cell_spike_times(cell_spike_times >= pre_laser_start & cell_spike_times <= pre_laser_stop);
        laser_spikes{this_cell,pulse} = cell_spike_times(cell_spike_times >= laser_start & cell_spike_times <= laser_stop);
        post_laser_spikes{this_cell,pulse} = cell_spike_times(cell_spike_times >= post_laser_start & cell_spike_times <= post_laser_stop);
        
        num_spikes_no_laser(this_cell,pulse) = numel(pre_laser_spikes{this_cell,pulse}) + numel(post_laser_spikes{this_cell,pulse});
        num_spikes_laser(this_cell,pulse) = numel(laser_spikes{this_cell,pulse});
        
        all_near_laser_spikes{this_cell,pulse} = [rel_pre_laser_spikes{this_cell,pulse}; rel_laser_spikes{this_cell,pulse}; rel_post_laser_spikes{this_cell,pulse}];
        temp_hist = hist(all_near_laser_spikes{this_cell,pulse},bins);
        temp_hist(temp_hist > 1) = 1;
        temp_hist(temp_hist == 0) = nan;
        all_near_laser_spikes_hist{this_cell,pulse} = temp_hist;
    end
end

firing_rate_laser_off = mean(num_spikes_no_laser,2)./(PRE_LASER_TIME/10^6 + POST_LASER_TIME/10^6);
firing_rate_laser_on = mean(num_spikes_laser,2)./(LASER_PULSE_TIME/10^6);
modulation_index = (firing_rate_laser_off - firing_rate_laser_on)./(firing_rate_laser_off + firing_rate_laser_on);

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

all_near_laser_spikes_hist_sorted_by_depth = all_near_laser_spikes_hist(sorted_depth_idx,:);
sorted_depth = depth(sorted_depth_idx);
sorted_modulation_value = modulation_value(sorted_depth_idx);

figure;
% cmap = jet(numCells);
cmap = flipud(jet(length(modulation_bins)-1));

for this_cell = 1:numCells
    for pulse = 1:numLaserPulses
%         y_value = numLaserPulses*(numChannels - sorted_depth(this_cell)) + pulse;
        y_value = numLaserPulses*(numCells + 1 - this_cell) + pulse - numLaserPulses;
        scatter(bins,all_near_laser_spikes_hist_sorted_by_depth{this_cell,pulse}*y_value,4,cmap(sorted_modulation_value(this_cell),:),'filled')
        hold on
    end
   
   
end
xlim([min(bins) max(bins)]) 
ylim([0 numLaserPulses*numCells])
title('Inhibition of MEC by Jaws')
xticks(-2.5:0.5:7.5)
xlabel('time, relative to laser onset (s)')
yticks(numLaserPulses/2:numLaserPulses:numLaserPulses*numCells-numLaserPulses/2)
yticklabels(numCells:-1:1)
ylabel('cell number')
axis tight

colormap(jet)
cb = colorbar;
set(cb,'YDir','reverse')
ylabel(cb,'Modulation Index')
caxis([-1 1])

yyaxis right
axr = gca;
axr.YLim = [0 numLaserPulses*numCells];
axr.YTick = numLaserPulses/2:numLaserPulses:numLaserPulses*numCells-numLaserPulses/2;
axr.YTickLabels = flipud(sorted_depth);
ylabel('Position Along Probe')
saveas(gcf,fullfile(kilosort_dir,'Raster Plot'),'png')

figure()
hist(modulation_index)
xlim([-1 1])
title('Modulation Index')
ylabel('Count')
xlabel('(FR Laser ON - FR Laser OFF)/(FR Laser ON + FR Laser OFF)')
saveas(gcf,fullfile(kilosort_dir,'Modulation Index'),'png')
