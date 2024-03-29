function filename = neuralynx2kilosortSherlock(directoryPath, target_dir, timeStamps)


    % convert H3 probe data to bin format for kilosort analysis
    % created by Kei Masuda
    % IL edited 4/9/19
    %
    % Inputs:
    % -------
    % directory : string
    %   contains path to data folder, no trailing backslash
    % target_dir : string
    %   path to folder to save output; optional, default is directorypath
    %
    % Outputs:
    % --------
    % .bin file that is a data matrix with samples of size 64channels x numOfSamples
    % if multiple files from the same date, they will be concatenated in
    %   order of recording
    
    % check input arguments
    if iscell(directoryPath)
        directoryPath = string(directoryPath);
    end
    if nargin == 1
        target_dir = directoryPath;
    end
    
    fprintf(strcat('\nStart Processing: ', datestr(now,'mmmm dd, yyyy HH:MM:SS AM'),'\n'));
    numOfChannels = 64;
%     timeStamps = '620031416955_620966217548';

    for csc = 1:numOfChannels
       
        keyboard
        cscPath = fullfile(directoryPath, ['CSC' num2str(csc) '_' timeStamps '.ncs']);

        % load neuralynx file, linearize samples, convert to int16
    %    [Samples,header]=Nlx2MatCSC(cscPath, [0 0 0 0 1], 1, 1, [] );
         [Samples,header]=Nlx2MatCSC_v3(cscPath, [0 0 0 0 1], 1, 1, [] );
        tmp=split(header{17}); %assuming conversion factor is in here;
        conv_factor = str2double(tmp{2}); % in volts
        conv_factor = conv_factor*10e6; % in micro volts
        %conv_factor = 1;
        Samples = reshape(Samples,1,[])*-1*conv_factor; 
        Samples= int16(Samples); 

        % pre-allocate dataMatrix size on first pass of each session
        if csc == 1 
           sample_idx = size(Samples,2);
           dataMatrix = zeros(numOfChannels,sample_idx,'int16');
        end

        % allocate csc data to appropriate rows/columns of data matrix
        dataMatrix(csc, :) = Samples; 
        fprintf(strcat('\nProcessed: ', num2str(csc), ' out of 64 CSC files.'));
    end    
    
    [~,Name,~] = fileparts(directoryPath);
    filename = [Name '_dataMatrix.bin'];
    fid = fopen(fullfile(target_dir, filename), 'w'); 
    fwrite(fid, dataMatrix, 'int16');
    fclose(fid);
    fprintf(strcat('\nDone Processing: ', datestr(now,'mmmm dd, yyyy HH:MM:SS AM'),'\n'));
end