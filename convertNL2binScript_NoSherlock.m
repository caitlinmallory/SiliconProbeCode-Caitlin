addpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin')
addpath('C:\Users\Caitlin Mallory\Documents\GitHub\SiliconProbeCode-Caitlin\MatlabImportExport_v6.0.0')
%%  

sourcedir = 'F:\2020-05-16_08-07-57\16352173949_17269640550';
timestamps = strsplit(sourcedir,'\'); timestamps = timestamps{end};

neuralynx2kilosort(sourcedir, sourcedir, timestamps)