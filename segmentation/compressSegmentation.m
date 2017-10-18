function compressSegmentation(p)
    % compressSegmentation(p)
    %   Compresses a segmentation, if possible.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % configuration
    resolutions = 2 .^ (0:9);
    taskCount = 50;
    
    % only WKW datasets can be compressed
    if ~isfield(p.seg, 'backend'); return; end
    if p.seg.backend ~= 'wkwrap'; return; end
    
    % find input root
    inRoot = buildInRoot(p.seg);
    bakRoot = strcat(inRoot, '-uncompressed');
    outRoot = strcat(inRoot, '-wip');
    
    % create output directory, if needed
    assert(not(exist(bakRoot, 'dir')));
    assert(not(exist(outRoot, 'dir')));
    mkdir(outRoot);
    
    for curRes = resolutions
        curInRoot = fullfile(inRoot, num2str(curRes));
        curOutRoot = fullfile(outRoot, num2str(curRes));
        
        if ~exist(curInRoot, 'dir')
            warning('Resolution %d not found → skipped');
            continue;
        end
        
        % prepare compressed WKW dataset
        wkwInit('compress', curInRoot, curOutRoot);
        
        fprintf('Compressing resolution %d ... ', curRes);
        wkwCompressDir(curInRoot, curOutRoot, taskCount);
        fprintf('✔\n');
    end
    
    % replace segmentation
    assert(movefile(inRoot, bakRoot));
    assert(movefile(outRoot, inRoot));
end

function inRoot = buildInRoot(inParam)
    dirParts = strsplit(inParam.root, filesep);
    
    % drop trailing file separator, if any
    if isempty(dirParts{end}); dirParts(end) = []; end
    
    % go one level up
    dirParts(end) = [];
    inRoot = strjoin(dirParts, filesep);
end
