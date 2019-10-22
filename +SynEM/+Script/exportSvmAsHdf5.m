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
    infoToHdf5(curOutFile, info);
end
