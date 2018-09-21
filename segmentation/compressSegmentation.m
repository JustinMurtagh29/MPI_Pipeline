function compressSegmentation(seg)
    % compressSegmentation(seg)
    %   Compresses a segmentation, if possible.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % configuration
    taskCount = 50;
    
    % only WKW datasets can be compressed
    if ~isfield(seg, 'backend'); return; end
    if seg.backend ~= 'wkwrap'; return; end
    
    % find input root
    inRoot = buildInRoot(seg);
    bakRoot = strcat(inRoot, '-uncompressed');
    outRoot = strcat(inRoot, '-wip');
    
    % create output directory, if needed
    assert(not(exist(bakRoot, 'dir')), ...
        'Folder %s already exists.', bakRoot);
    assert(not(exist(outRoot, 'dir')), ...
        'Folder %s already exists.', outRoot);
    mkdir(outRoot);
    
    resDirs = dir(inRoot);
    resDirs = {resDirs([resDirs.isdir]).name};
    
    % Only keep directories with valid names:
    % That is, either {:num:} or {:num:}-{:num:}-{:num:}
    resDirs(cellfun(@isempty, regexp( ...
        resDirs, '^\d+(-\d+-\d+)?$'))) = [];
    
    for curIdx = 1:numel(resDirs)
        curRes = resDirs{curIdx};
        
        curInRoot = fullfile(inRoot, curRes);
        curOutRoot = fullfile(outRoot, curRes);
        
        % prepare compressed WKW dataset
        wkwInit('compress', curInRoot, curOutRoot);
        
        fprintf('Compressing resolution %s ... ', curRes);
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
