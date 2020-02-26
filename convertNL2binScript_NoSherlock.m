addpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin')
addpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin\MatlabImportExport_v6.0.0')
%%  

sourcedir = 'F:\2020-02-24_13-25-55\19566053356_20051634130';
timestamps = strsplit(sourcedir,'\'); timestamps = timestamps{end};

neuralynx2kilosort(sourcedir, sourcedir, timestamps)