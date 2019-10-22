% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
inDir = '/gaba/u/bstaffle/data/CNN_Training/SVM_TrainingCubes/';
outDir = '/home/amotta/Desktop';

cubeIds = 1:7;
segCatVals = uint32(0:3);
segCatNames = {'Background', 'Synapse', 'Vesicle', 'Mitochondrion'};

info = Util.runInfo();
Util.showRunInfo(info);

%% Build info string
clear cur*;
curIsDirtyStrs = {'', ' (dirty)'};

infoStr = cell(0, 1);
infoStr{end + 1} = sprintf('%s', info.filename);

for curRepoId = 1:numel(info.git_repos)
    curRepo = info.git_repos{curRepoId};
    infoStr{end + 1} = sprintf( ...
        '%s %s%s', curRepo.remote, curRepo.hash, ...
        curIsDirtyStrs{2 - isempty(curRepo.diff)}); %#ok
end

infoStr{end + 1} = sprintf( ...
    '%s@%s. MATLAB %s. %s', ...
    info.user, info.hostname, ...
    info.matlab_version, info.time);
infoStr = strjoin(infoStr, newline);

%% Export cubes
clear cur*;

for curId = cubeIds
    curInFile = fullfile(inDir, sprintf('SVM_%d.mat', curId));
    curOutFile = fullfile(outDir, sprintf('cube-%d.hdf5', curId));
    
    curIn = load(curInFile);
    curIn.seg = categorical(curIn.seg, segCatVals, segCatNames);
    
    numericToHdf5(curOutFile, '/em', curIn.raw);
    categoricalToHdf5(curOutFile, '/label', curIn.seg);
    h5writeatt(curOutFile, '/', 'info', infoStr);
end

%% Utilities
function categoricalToHdf5(outFile, dset, cats)
    assert(iscategorical(cats));
    catNames = categories(cats);
    
    assert(numel(catNames) <= intmax('uint8'));
    numericToHdf5(outFile, dset, uint8(cats));
    
    for curIdx = 1:numel(catNames)
        curName = lower(catNames{curIdx});
        h5writeatt(outFile, dset, curName, uint8(curIdx));
    end
end

function numericToHdf5(file, dset, data)
    assert(isnumeric(data));
    sz = size(data);
    
    % NOTE(amotta): Remove trailing singleton dimension
    if numel(sz) == 2 && sz(2) == 1; sz = sz(1); end
    
    % NOTE(amotta): Compression is counter-productive in this particular
    % case. Storage usage reported by h5ls is around 30 % to 50 % with
    % deflate (9).
    h5create(file, dset, sz, 'Datatype', class(data));
    h5write(file, dset, data);
end
