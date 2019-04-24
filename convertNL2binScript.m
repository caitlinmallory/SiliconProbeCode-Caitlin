addpath('/home/users/cmallory/SiliconProbeCode-Caitlin')
addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/MatlabImportExport_v6.0.0')
addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/releaseDec2015')

%%  
sourcedir ='/oak/RightHemProbeRecording';
targetdir = '/oak/RightHemProbeRecording/Processed_Data';


neuralynx2kilosortSherlock(sourcedir, targetdir)
