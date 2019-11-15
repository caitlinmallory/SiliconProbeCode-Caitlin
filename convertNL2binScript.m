addpath('/home/users/cmallory/SiliconProbeCode-Caitlin')
%addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/MatlabImportExport_v6.0.0')
addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/releaseDec2015')
addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/releaseDec2015/binaries')
addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/releaseDec2015/source')
%%  
sourcedir ='/home/users/cmallory/oak/11_11_19/18794456157_18990710047';
targetdir = '/home/users/cmallory/oak/RightHemProbeRecording/11_11_19/18794456157_18990710047/Processed_Data';


neuralynx2kilosortSherlock(sourcedir, targetdir)
