% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_classified.mat');
outDir = '/home/amotta/Desktop';

% List of target class to export
% One NML file per target class will be generated
exportClasses = {'Smoothies', 'Apicals'};
exportMaxTrees = 500;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

dend = load(dendFile);

%% Generating NML files
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curClassIdx = 1:numel(exportClasses)
    curClass = exportClasses{curClassIdx};
    
    curDendIds = dend.(sprintf('ind%s', curClass));
    curDendIds = find(curDendIds & dend.indBigDends);
    
    % Sort by size
    curSortIds = arrayfun( ...
        @(dend) size(dend.nodes, 1), ...
        dend.dendrites(curDendIds));
   [~, curSortIds] = sort(curSortIds, 'descend');
    curDendIds = curDendIds(curSortIds);
    
    % Limit
    curTreeCount = min(numel(curDendIds), exportMaxTrees);
    curDendIds = curDendIds(1:curTreeCount);
    
    curSkel = skel;
    for curId = reshape(curDendIds, 1, [])
        curAgglo = dend.dendrites(curId);
        curSkel = curSkel.addTree( ...
            sprintf('Dendrite %d', curId), ...
            curAgglo.nodes(:, 1:3), curAgglo.edges);
    end
    
    curSkelName = sprintf('%s.nml', lower(curClass));
    curSkel.write(fullfile(outDir, curSkelName));
end
