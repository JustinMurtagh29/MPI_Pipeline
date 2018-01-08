% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%%
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
thisDir = fileparts(mfilename('fullpath'));

nmlDir = fullfile( ...
    thisDir, '..', '..', ...
    'evaluationData', 'new_axon_gt_ROI2017');

%% load parameters
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

typeProbs = struct;
typeProbs.maxSegId = Seg.Global.getMaxSegId(param);
typeProbs = connectEM.addSegmentClassInformation(param, typeProbs);

%% find NML files
nmlFiles = dir(fullfile(nmlDir, '*.nml'));
nmlFiles = fullfile(nmlDir, reshape({nmlFiles.name}, [], 1));

%% load
nmlNodes = cell(size(nmlFiles));

for curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    
    curNml = slurpNml(curNmlFile);
    curNodes = NML.buildNodeTable(curNml);
    
    curNodes.coord = curNodes.coord + 1;
    curNodes.segId(:) = -1;
    
    curMask = ...
        all(curNodes.coord >= param.bbox(:, 1)', 2) ...
      & all(curNodes.coord <= param.bbox(:, 2)', 2);
    curNodes.segId(curMask) = ...
        Seg.Global.getSegIds(param, curNodes.coord(curMask, :));
    
    nmlNodes{curIdx} = curNodes;
end

%% show results
figure;

for curIdx = 1:numel(nmlFiles)
    curSegIds = nmlNodes{curIdx}.segId;
    curSegIds = curSegIds(curSegIds > 0);
    curAxonProbs = typeProbs.axonProb(curSegIds);
    
    subplot(numel(nmlFiles), 1, curIdx);
    histogram(curAxonProbs, linspace(0, 1, 21));
end