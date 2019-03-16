% NOTE: This script was written in such a way that additional query rounds
% can easily be implemented by replacing the paths in axonNmlFiles with the
% results from the previous annotation round. Slight changes to the query
% selection section should be all that's neeeded to get this going.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
% NOTE(amotta): This file is identical to 20180726T190355_results.mat with
% the exception of the additional `outputMap.axonData.segIds` field.
outputConnFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';
outputConnRunId = '20190221T112510';

annDir = 'plass-proofreading_round-1';
annDir = connectEM.WholeCell.Data.getFile(annDir);

queryCount = 0;
queryDir = '/home/amotta/Desktop';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

% Loading tracings and synapses of L4 cells
[curDir, curFile] = fileparts(outputConnFile);
curFile = sprintf('%s__%s_connectome.mat', curFile, outputConnRunId);
curFile = fullfile(curDir, curFile);

curData = load(curFile);
synT = curData.synT;

outputMapFile = curData.info.param.outputMapFile;
curData = load(outputMapFile);
axonData = curData.axonData;

connFile = curData.info.param.connFile;
[~, synAgglos] = connectEM.Connectome.load(param, connFile);

synAgglos = cellfun(@union, ...
    synAgglos.synapses.presynId, ...
    synAgglos.synapses.postsynId, ...
    'UniformOutput', false);

%% Process existing annotations
clear cur*;

annT = {'axonId', 'synId', 'isCorrect', 'isSpine', 'somaDist'};
annT = array2table(nan(0, 5), 'VariableNames', annT);
annT.isCorrect = logical(annT.isCorrect);
annT.isSpine = logical(annT.isSpine);

if ~isempty(annDir)
    for curAxonId = 1:numel(axonData)
        curAxonData = axonData(curAxonId);
        curSyn = curAxonData.synapses;
        
       [~, curNmlFile] = fileparts(curAxonData.nmlFile);
        curNmlFile = fullfile(annDir, sprintf('%s.nml', curNmlFile));
        curAnn = skeleton(curNmlFile);
        
        curAnn = curAnn.names;
        curAnn = curAnn(startsWith(curAnn, 'Synapse'));
        if isempty(curAnn); continue; end
        
        curIsCorrect = contains(curAnn, 'correct', 'IgnoreCase', true);
        curSynIds = regexp(curAnn, '^Synapse (\d+)', 'tokens', 'once');
        assert(not(any(cellfun(@isempty, curSynIds))));
        
        curSynIds = cat(1, curSynIds{:});
        curSynIds = cellfun(@str2double, curSynIds);
        
       [~, curIsSpine] = ismember(curSynIds, synT.id);
        curIsSpine = synT.isSpine(curIsSpine);
        
       [~, curSomaDist] = ismember(curSynIds, curSyn.id);
        curSomaDist = curSyn.somaDist(curSomaDist);
        
        curAnnT = table;
        curAnnT.synId = curSynIds;
        curAnnT.isCorrect = curIsCorrect;
        curAnnT.isSpine = curIsSpine;
        curAnnT.somaDist = curSomaDist;
        curAnnT.axonId(:) = curAxonId;
        
        annT = cat(1, annT, curAnnT);
    end
end

%% Visualize annotations
clear cur*;

% Quantitative
curEvalT = table;
curEvalT.synType = {'Spine'; 'Shaft'};
curEvalT.somaDistsUm = { ...
    annT.somaDist(annT.isCorrect & annT.isSpine) / 1E3; ...
    annT.somaDist(annT.isCorrect & not(annT.isSpine)) / 1E3};

curEvalT.mean = cellfun(@mean, curEvalT.somaDistsUm);
curEvalT.std = cellfun(@std, curEvalT.somaDistsUm);

fprintf('PLASS evaluation\n\n');
disp(curEvalT); fprintf('\n');

[~, curPVal] = ttest2(curEvalT.somaDistsUm{:});
fprintf('* t-test: p = %f\n', curPVal);

curPVal = ranksum(curEvalT.somaDistsUm{:});
fprintf('* rank sum test: p = %f\n', curPVal);

% Visual
curBinEdges = linspace(0, 200, 11);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

histogram( ...
    curAx, curEvalT.somaDistsUm{1}, ...
    'BinEdges', curBinEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'EdgeColor', 'magenta');
histogram( ...
    curAx, curEvalT.somaDistsUm{2}, ...
    'BinEdges', curBinEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'EdgeColor', 'black');

curLeg = {'Spine synapses', 'Shaft synapses'};
curLeg = legend(curAx, curLeg, 'Location', 'SouthOutside');

xticks(curAx, curBinEdges(1):10:curBinEdges(end));
curAx.XTickLabel(mod(0:numel(curAx.XTickLabel), 5) ~= 0) = {''};
xlabel(curAx, 'Axonal path length from soma [µm]');
ylabel(curAx, 'Probability');

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [250, 235];

%% Select synapses to query
clear cur*;
rng(0);

queries = struct;
queries.spine = find(synT.isSpine);
queries.shaft = find(not(synT.isSpine));

% Shuffle
queries = structfun( ...
    @(v) v(randperm(numel(v))), ...
    queries, 'UniformOutput', false);

% Select subset
queries = structfun( ...
    @(v) v(1:min(numel(v), queryCount)), ...
    queries, 'UniformOutput', false);

%% Generate NML files
clear cur*;
curSynClasses = {'Shaft', 'Spine'};

for curAxonId = 1:numel(axonData)
    curAxonNmlFile = axonData(curAxonId).nmlFile;
   [~, curAxonName] = fileparts(curAxonNmlFile);
    
    curSynT = queries;
    curSynT = structfun( ...
        @(v) v(synT.preAggloId(v) == curAxonId), ...
        curSynT, 'UniformOutput', false);
    
    curSynT = cell2mat(struct2cell(curSynT));
    curSynT = synT(curSynT, :);
    
    curSkel = skeleton(curAxonNmlFile);
    curSkel = curSkel.deleteTreeWithName('Dendrite');
    
    for curSynIdx = 1:height(curSynT)
        curSyn = curSynT(curSynIdx, :);
        curSegIds = synAgglos{curSyn.id};
        
        curSkel = Skeleton.fromMST( ...
            segPoints(curSegIds, :), ...
            param.raw.voxelSize, curSkel);
        curSkel.names{end} = sprintf( ...
            'Synapse %d. %s. Todo', curSyn.id, ...
            curSynClasses{1 + curSyn.isSpine});
    end
    
    curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
    curSkel.write(fullfile(queryDir, sprintf('%s.nml', curAxonName)));
end
