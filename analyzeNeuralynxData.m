addpath(genpath('MatlabImportExport_v6.0.0'))
Timestamps = cell(64,1);
Samples = cell(64,1);
Features = cell(65,1);
data_path = 'F:\Foster Lab Data\2019-08-09_12-11-00\6159665707-7893828980\';
ADBitVolts = 0.000000015258789062500000;
recording_times = '_6159665707_7893828980';


for c = 1:64
    
[Timestamps{c}, ~, ~, Features{c}, Samples{c}, Header] = Nlx2MatSpike(strcat(data_path,'SE',num2str(c), recording_times, '.nse'), [1 1 1 1 1], 1, 1, [] );

end

% for each electrode, find when then cell crossed a threshold
spike_threshold = 50;
Good_spikes = cell(64,1);

for c = 1:64
    
    good_spikes = [];
    cell_times = Timestamps{c};
    
    for t = 1:length(Samples{c})
            max_value = max(Samples{c}(:,:,t).*ADBitVolts*10^6);
            if max_value > spike_threshold
            good_spikes = [good_spikes; cell_times(t)];
            end
    end
     Good_spikes{c} = good_spikes;   
end
    
event_file = fullfile(data_path,'Events_6159665707_7893828980.csv');
events = readtable(event_file);


laser_start = events.TimeStamp(events.TTLValue == 1);
laser_duration = 10^6; %1 sec
laser_stop = laser_start + laser_duration;

laser_off_start = laser_stop;
laser_off_stop = laser_off_start + 9*10^6;

spike_count_laser_on = nan(64,1);
spike_count_laser_off = nan(64,1);
spike_ind_laser_on = cell(64,1);
spike_ind_laser_off = cell(64,1);

for c = 1:64
    total_spike_count_laser_on = 0;
    all_spikes_laser_on = [];
    for x = 1:length(laser_start)
        new_spikes_laser_on = find(Good_spikes{c} >= laser_start(x) & Good_spikes{c} <= laser_stop(x));
        new_spike_count_laser_on = numel(find(Good_spikes{c} >= laser_start(x) & Good_spikes{c} <= laser_stop(x)));
        total_spike_count_laser_on = total_spike_count_laser_on + new_spike_count_laser_on;
        all_spikes_laser_on = [all_spikes_laser_on; new_spikes_laser_on];
    end
    spike_ind_laser_on{c} = all_spikes_laser_on;
    spike_count_laser_on(c) = total_spike_count_laser_on;
end

for c = 1:64
    total_spike_count_laser_off = 0;
    all_spikes_laser_off = [];
    for x = 1:length(laser_start)
        new_spikes_laser_off = find(Good_spikes{c} >= laser_off_start(x) & Good_spikes{c} <= laser_off_stop(x));
        new_spike_count_laser_off = numel(find(Good_spikes{c} >= laser_off_start(x) & Good_spikes{c} <= laser_off_stop(x)));
        total_spike_count_laser_off = total_spike_count_laser_off + new_spike_count_laser_off;
        all_spikes_laser_off = [all_spikes_laser_off; new_spikes_laser_off];
    end
    spike_ind_laser_off{c} = all_spikes_laser_off;
    spike_count_laser_off(c) = total_spike_count_laser_off;
end

fr_laser_on = spike_count_laser_on./(1*length(laser_start));
fr_laser_off = spike_count_laser_off./(9*length(no_laser_start));


plot(Samples{1,1}(:,:,68).*ADBitVolts*10^6)




%   Feature 1 = Peak
%   Feature 2 = Valley
%   Feature 3 = Energy
%   Feature 4 = Height
%   
%     {'-Feature Peak 0 0 0 31 1'                                      }
%     {'-Feature Valley 1 0 0 31 1'                                    }
%     {'-Feature Energy 2 0 0 31 1'                                    }
%     {'-Feature Height 3 0 0 31 1'                                    }
%     {'-Feature NthSample 4 0 0 31 1 4'                               }
%     {'-Feature NthSample 5 0 0 31 1 16'                              }
%     {'-Feature NthSample 6 0 0 31 1 24'                              }
%     {'-Feature NthSample 7 0 0 31 1 28'                              }