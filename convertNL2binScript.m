addpath('/home/users/cmallory/SiliconProbeCode-Caitlin')
% addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/MatlabImportExport_v6.0.0')
addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/releaseDec2015')
addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/releaseDec2015/binaries')
addpath('/home/users/cmallory/SiliconProbeCode-Caitlin/releaseDec2015/source')
%%  
sourcedir ='/home/users/cmallory/oak/11_18_19/620031416955_620365221790/';
targetdir ='/home/users/cmallory/oak/11_18_19/620031416955_620365221790/';
timestamps = '620031416955_620365221790';
% sourcedir = 'Z:\giocomo\cmallory\11_18_19\616169768506_616432328369';
% targetdir = 'Z:\giocomo\cmallory\11_18_19\616169768506_616432328369';

neuralynx2kilosortSherlock(sourcedir, targetdir, timestamps)