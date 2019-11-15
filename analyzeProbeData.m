% analyzeProbeData
% cids = cluster group; 1 = MUA, 2 = Good, 3 = Unsorted
% sp = spike times
% sp.clu = clusterID

addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\Spikes'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\npy-matlab'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin'))
addpath(genpath('C:\Users\Caitlin Mallory\Documents\GitHub\MatlabImportExport_v6.0.0'));
%myKsDir = 'F:\Foster Lab Data\2018-12-19_11-34-45\RightHemProbeRecording\preAutoMerge';
myKsDir = 'F:\Foster Lab Data\2019-08-09_12-11-00\6159665707-7893828980\kilosort';

% Neuralynx data path and session info:
myNlxDir = 'F:\Foster Lab Data\2019-08-09_12-11-00\6159665707-7893828980\';
Nlx_recording_times = '_6159665707_7893828980';
Nlx_lfp_channels = 1:64;

cfg_def.fc = {};
cfg_def.TimeConvFactor = 1; % 10^-6 means convert nlx units to seconds
cfg_def.VoltageConvFactor = 10^6; % 1 means output in volts, 1000 in mV, 10^6 in uV
%%
chunk_time = 3; %sec
Fs = 32000;
mean_lfp_data_near_laser = nan(length(Nlx_lfp_channels), Fs*chunk_time);
mean_lfp_data_near_laser_normalized = nan(length(Nlx_lfp_channels), Fs*chunk_time);

for c = 1:length(Nlx_lfp_channels)
% load a neuralynx eeg file to get the timestamps
[lfp_ts, ~, sample_frequencies, number_of_valid_samples, lfp_samples, header] = Nlx2MatCSC(strcat(myNlxDir,'CSC',num2str(Nlx_lfp_channels(c)),Nlx_recording_times,'.ncs'),[1 1 1 1 1],1,1,[]);

% check for constant sampling frequency
Fs = unique(sample_frequencies);
if length(Fs) ~= 1
    error('More than one sampling frequency found!');
end

% extract information from header
hdr = readCSCHeader(header);

% convert units
lfp_samples = lfp_samples .*cfg_def.VoltageConvFactor .*hdr.ADBitVolts;

% construct within-block tvec
nSamplesPerBlock = size(lfp_samples,1);
block_tvec = (0:(1/Fs):((nSamplesPerBlock-1)/Fs)).*10^6;

% slow but tractable: loop over blocks remembering to initialize variables
% first

lfp_samples_linearized = nan(numel(lfp_samples),1); % allocate memory and then fill; can trim later
lfp_ts_linearized = nan(numel(lfp_samples),1);

idx = 1; % move this along as we go through the loop over all samples

nBlocks = length(lfp_ts);
bad_blocks = 0; % counter
for iB = 1:nBlocks
    
    nvs = number_of_valid_samples(iB);
    if nvs ~= 512, bad_blocks = bad_blocks + 1; end
    
    current_data = lfp_samples(1:nvs,iB);
    current_time = block_tvec(1:nvs)+lfp_ts(iB);
    
    lfp_samples_linearized(idx:idx+nvs-1) = current_data;
    lfp_ts_linearized(idx:idx+nvs-1) = current_time;
    
    idx = idx + nvs;
    
end % of block loop

cfg.badBlocks = badBlocks;

% load the events file
event_file = fullfile(myNlxDir,strcat('Events',Nlx_recording_times,'.csv'));
events = readtable(event_file);
session_length = (max(lfp_ts) - min(lfp_ts))/10^6; %sec

laser_1000mW_start = find(strcmp(events.Eventstring,'Increasing laser power to 1000mA (Full power)'));
laser_1000mW_stop = find(strcmp(events.Eventstring,'Lowering laser power to 750mA'));

events_1000mW = events(laser_1000mW_start:laser_1000mW_stop,:);

laser_start = events_1000mW.TimeStamp(events_1000mW.TTLValue == 1);
laser_off_start  = events_1000mW.TimeStamp(events_1000mW.TTLValue == 0);
laser_off_start(1) = [];
laser_off_start(end) = [];
laser_stop = laser_off_start - (1/Fs)*10^6;
laser_off_stop = laser_start(2:end) - (1/Fs)*10^6;
laser_off_stop = [laser_off_stop; events_1000mW.TimeStamp(strcmp(events_1000mW.Eventstring,'Lowering laser power to 750mA'))];
pre_laser_time = 10^6; %us
laser_time = 10^6; %us
post_laser_time = 10^6; %us

lfp_data_near_laser = nan(length(laser_start),Fs*chunk_time);
lfp_data_near_laser_normalized = nan(length(laser_start),Fs*chunk_time);

for x = 1:length(laser_start)
     [~,ind_laser_start] = min(abs(lfp_ts_linearized - laser_start(x)));
     ind_seg_start = ind_laser_start - pre_laser_time*Fs/10^6;
     
     [~,ind_laser_stop] = min(abs(lfp_ts_linearized - laser_stop(x)));
     ind_seg_stop = ind_laser_stop + post_laser_time*Fs/10^6;
    
     lfp_samples_linearized_laser = lfp_samples_linearized(ind_laser_start:ind_laser_stop);
     
     % this is supposed to be 32000 samples (1 sec) but the actual laser ON times are inconsisent.
     % Interpolate between the start and stop of the laser.
     
     x_samples = 1:numel(lfp_samples_linearized_laser);
     xq_samples = linspace(1,numel(lfp_samples_linearized_laser),32000);
     lfp_samples_linearized_laser_interp = interp1(x_samples,lfp_samples_linearized_laser,xq_samples);
     
     lfp_data_near_laser(x,:) = [lfp_samples_linearized(ind_seg_start:ind_laser_start-1)' lfp_samples_linearized_laser_interp lfp_samples_linearized(ind_laser_stop + 1: ind_seg_stop)'];
     lfp_data_near_laser_normalized(x,:) = (lfp_data_near_laser(x,:)-(min(lfp_data_near_laser(x,:))))/(max(lfp_data_near_laser(x,:)-min(lfp_data_near_laser(x,:))));
end

mean_lfp_data_near_laser(c,:) = mean(lfp_data_near_laser,1);
mean_lfp_data_near_laser_normalized(c,:) = mean(lfp_data_near_laser_normalized,1);
end
%%
channels = 1:9:65;
for c = 1:length(channels)
subplot(length(channels),1,c)
hold on
plot([size(mean_lfp_data_near_laser_normalized,2)/3 size(mean_lfp_data_near_laser_normalized,2)/3],[min(min(mean_lfp_data_near_laser_normalized(channels,:), [], 2)) max(max(mean_lfp_data_near_laser_normalized(channels,:), [], 2))],'-k')
hold on
plot([size(mean_lfp_data_near_laser_normalized,2)/3*2 size(mean_lfp_data_near_laser_normalized,2)/3*2],[min(min(mean_lfp_data_near_laser_normalized(channels,:), [], 2)) max(max(mean_lfp_data_near_laser_normalized(channels,:), [], 2))],'-k')
hold on
plot(mean_lfp_data_near_laser_normalized(channels(c),:),'r','LineWidth',2)
ylabel(['channel ' num2str(channels(c))])
end




% load the silicon probe cluster_info file containing basic info about each
% clusete
% remove noise clusters
[~,header,raw] = tsvread( fullfile(myKsDir,'cluster_info.tsv'));
noise_clusters = find(strcmp(raw(2:end,6),'noise'));
depth = raw(2:end,4);
depth = cellfun(@str2num, depth);
depth(noise_clusters) = [];

% load the phy output containing timestamps for each cluster
sp = loadKSdir(myKsDir);
spike_times = sp.st.*10^6; %change to microsec to be consistent with neuralynx
spike_times = spike_times + lfp_ts_linearized(1);

% pull out the cells
cluster_id = sp.cids;
good_cluster_id = cluster_id(sp.cgs==2);

% find the depth of the good cells
good_cluster_depth = depth(sp.cgs==2);


numCells = numel(good_cluster_id);

spike_count_laser_on = nan(numCells,1);
spike_count_laser_off = nan(numCells,1);
spike_ind_laser_on = cell(numCells,1);
spike_ind_laser_off = cell(numCells,1);

sessionFR = nan(numCells,1);
restrictedFR = nan(numCells,1);
rangeSpiking = nan(numCells,1);

for k = 1:numCells
    restrictedFR(k) = numel(sp.st(sp.clu == good_cluster_id(k)))/(max(sp.st(sp.clu == good_cluster_id(k)))-min(sp.st(sp.clu == good_cluster_id(k))));
    sessionFR(k) = numel(sp.st(sp.clu == good_cluster_id(k)))/session_length;
    firstSpike =  min(sp.st(sp.clu == good_cluster_id(k)));
    lastSpike = max(sp.st(sp.clu == good_cluster_id(k)));
    rangeSpiking(k) = lastSpike-firstSpike;
end

subplot(1,2,1)
boxplot(sessionFR)
ylabel('spikes/sec');
title('Session FR')

subplot(1,2,2)
boxplot(rangeSpiking)
ylabel('sec');
title(['time active; session = ',num2str(session_length), 'sec'])


for c = 1:numCells
    total_spike_count_laser_on = 0;
    all_spikes_laser_on = [];
    for x = 1:length(laser_start)
        cell_spike_times = spike_times(sp.clu == good_cluster_id(c));
        new_spikes_laser_on = cell_spike_times(cell_spike_times  >= laser_start(x) & cell_spike_times  <= laser_stop(x));
        new_spike_count_laser_on = numel(new_spikes_laser_on);
        total_spike_count_laser_on = total_spike_count_laser_on + new_spike_count_laser_on;
        all_spikes_laser_on = [all_spikes_laser_on; new_spikes_laser_on];
    end
    spike_ind_laser_on{c} = all_spikes_laser_on;
    spike_count_laser_on(c) = total_spike_count_laser_on;
end

for c = 1:numCells
    total_spike_count_laser_off = 0;
    all_spikes_laser_off = [];
    for x = 1:length(laser_off_start)
        cell_spike_times = spike_times(sp.clu == good_cluster_id(c));
        new_spikes_laser_off = cell_spike_times(cell_spike_times  >= laser_off_start(x) & cell_spike_times  <= laser_off_stop(x));
        new_spike_count_laser_off = numel(new_spikes_laser_off);
        total_spike_count_laser_off = total_spike_count_laser_off + new_spike_count_laser_off;
        all_spikes_laser_off = [all_spikes_laser_off; new_spikes_laser_off];
    end
    spike_ind_laser_off{c} = all_spikes_laser_off;
    spike_count_laser_off(c) = total_spike_count_laser_off;
end

fr_laser_on = spike_count_laser_on./(1*length(laser_start));
fr_laser_off = spike_count_laser_off./(9*length(laser_off_start));
fr_diff_normalized = (fr_laser_on-fr_laser_off)./fr_laser_off;

colorstring = 'bgrcmykbgrcm';
for c = 1:numCells
    plot([1 2],[fr_laser_off(c) fr_laser_on(c)],'-o','Color',colorstring(c),'MarkerFaceColor',colorstring(c))
    hold on
end
xlim([0 3])
ylabel('firing rate, Hz','fontSize',20)
xticks([1 2])
xticklabels({'Laser Off', 'Laser ON'})
