addpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin')
addpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin\MatlabImportExport_v6.0.0')
%%  

sourcedir = 'F:\2020-03-02_12-13-08\7486602027_7993394792';
timestamps = strsplit(sourcedir,'\'); timestamps = timestamps{end};

neuralynx2kilosort(sourcedir, sourcedir, timestamps)