function csc_info = readCSCHeader(Header)

csc_info = [];
for hline = 1:length(Header)
    
    line = strtrim(Header{hline});
    
    if isempty(line) || ~strcmp(line(1),'-') % not an informative line, skip
        continue;
    end
    
    % if we are reading an older version header with colons, remove them
    colon_idx = strfind(line,':');
    if ~isempty(colon_idx)
       line(colon_idx) = []; 
    end
    
    % This expression captures the first chunk of alphanumeric characters and puts it into
    % <key> then puts whatever is to the right of it into <val>. If there is only one
    % character chunk (e.g., missing value) then it returns <val> as empty.
    a = regexp(line(2:end),'(?<key>^\S+)\s+(?<val>.*)|(?<key>\S+)','names');

    % deal with characters not allowed by MATLAB struct
    if strcmp(a.key,'DspFilterDelay_ï¿½s') || strcmp(a.key,'DspFilterDelay_µs')
        a.key = 'DspFilterDelay_us';
    end
    
    csc_info = setfield(csc_info,a.key,a.val);
    
    % convert to double if possible
    if ~isnan(str2double(a.val))
        csc_info = setfield(csc_info,a.key,str2double(a.val));
    end
    
end