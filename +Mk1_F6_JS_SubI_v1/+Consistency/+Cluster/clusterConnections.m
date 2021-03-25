% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by 
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat');

info = Util.runInfo();
Util.showRunInfo(info);

asiRunId = '20210325T114245'; % NOTE: Using proofread sasd pairs

info = Util.runInfo();
Util.showRunInfo(info);

%% Load axon-spine interfaces
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = sprintf('%s__%s_sasdData.mat', curAsiFile, asiRunId);
curAsiFile = fullfile(curDir, curAsiFile);

m = load(curAsiFile);
dataTable = m.dataTable;

% choose which synapses to analyse
sasdType = 'spine-spine'; %

switch sasdType
    case 'spine-spine'
        synapseSearchList1 = {'Spine','Prim','Second'};
        synapseSearchList2 = {'Spine','Prim','Second'};
    case 'shaft-shaft'
        synapseSearchList1 = {'Shaft'};
        synapseSearchList2 = {'Shaft'};
    case 'spine-shaft'
        synapseSearchList1 = {'Spine','Prim','Second'};
        synapseSearchList2 = {'Shaft'};
    otherwise
        error('Pls specify correct sasd synapse types')
end

% SASD pairs
idxPlot = (contains(dataTable.syn1,synapseSearchList1) & contains(dataTable.syn2,synapseSearchList2)) | ...
          (contains(dataTable.syn1,synapseSearchList2) & contains(dataTable.syn2,synapseSearchList1)) ;
x1 = dataTable(idxPlot,:).asi1;
x2 = dataTable(idxPlot,:).asi2;

d = sort([x1,x2],2,'descend'); % sort asi data
x1 = d(:,1); x2 = d(:,2);
clear d idxPlot
assert(numel(x1) == numel(x2))
sprintf('Found %d SASD %s pairs',numel(x1), sasdType)

%{
asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

%% Prepare data
clear cur*;

plotConfigs = struct('synIds', {}, 'title', {}, 'tag', {});

axonClasses = { ...
    'Exc', {'Corticocortical', 'Thalamocortical'}; ...
    'CC',  {'Corticocortical'}; ...
    'TC',  {'Thalamocortical'}};
targetClasses = { ...
    'All', categories(asiT.targetClass); ...
    'PD',  'ProximalDendrite'; ...
    'AD',  'ApicalDendrite'; ...
    'OD',  'OtherDendrite'};

for curAxonIdx = 1:size(axonClasses, 1)
    curAxonClass = axonClasses(curAxonIdx, :);
    
    for curTargetIdx = 1:size(targetClasses, 1)
        curTargetClass = targetClasses(curTargetIdx, :);
        
        curSynIds = find( ...
            asiT.type == 'PrimarySpine' ...
          & ismember(asiT.axonClass, curAxonClass{2}) ...
          & ismember(asiT.targetClass, curTargetClass{2}));
        if isempty(curSynIds); continue; end
        
        curTitle = sprintf( ...
            '%s → %s primary spine synapses', ...
            curAxonClass{1}, curTargetClass{1});
        curTag = sprintf('%s %s pri sp', ...
            curAxonClass{1}, curTargetClass{1});
        
        curPlotConfig = plotConfigs([]);
        curPlotConfig(1).synIds = curSynIds;
        curPlotConfig(1).title = curTitle;
        curPlotConfig(1).tag = curTag;
        
        plotConfigs(end + 1) = curPlotConfig; %#ok
    end
end

plotConfigs = reshape( ...
    plotConfigs, ...
    size(targetClasses, 1), ...
    size(axonClasses, 1));
%}

%% Run model
clear cur*;
%{
curSaSdSynIds = ...
    connectEM.Consistency.buildPairConfigs(asiT, plotConfigs(1));
curSaSdSynIds = curSaSdSynIds(1).synIdPairs;
%}

% SL: add already extracted sasd pairs areas
asiT = struct;
asiT.area = [x1, x2]; % pairs

curH5File = strcat(tempname, '.h5');
for curIdx = 1:2
    curDset = sprintf('/log10Asi%d', curIdx);
%    curData = log10(asiT.area(curSaSdSynIds(:, curIdx)));   
    curData  = log10(asiT.area(:,curIdx));

    h5create(curH5File, curDset, size(curData), 'Datatype', class(curData));
    h5write(curH5File, curDset, curData);
end

curCmd = fileparts(mfilename('fullpath'));
curCmd = fullfile(curCmd, 'clusterConnections.py');
curCmd = sprintf('%s "%s"', curCmd, curH5File);
[curExit, curResult] = system(curCmd, '-echo');
