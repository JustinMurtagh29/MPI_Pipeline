function compressSegmentation(inParam, outRoot)
    % compressSegmentation(inParam, outRoot)
    %   Compresses a segmentation, if possible.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % only WKW datasets can be compressed
    if ~isfield(inParam, 'backend'); return; end;
    if inParam.backend ~= 'wkwrap'; return; end;
    
    inRoot = buildInRoot(inParam);
    resolutions = 2 .^ (0:9);
    taskCount = 50;
    
    for curRes = resolutions
        curInRoot = fullfile(inRoot, num2str(curRes));
        curOutRoot = fullfile(outRoot, num2str(curRes));
        
        if ~exist(curInRoot, 'dir')
            warning('Resolution %d not found → skipped');
            continue;
        end
        
        % prepare compressed WKW dataset
        wkwInit('compress', curInRoot, curOutRoot);
        
        fprintf('Compressing resolution %d... ', curRes);
        wkwCompressDir(curInRoot, curOutRoot, taskCount);
        fprintf('✔\n');
    end
end

function inRoot = buildInRoot(inParam)
    dirParts = strsplit(inParam.root, filesep);
    
    % drop trailing file separator, if any
    if isempty(dirParts{end}); dirParts(end) = []; end;
    
    % go one level up
    dirParts(end) = [];
    inRoot = strjoin(dirParts, filesep);
end
