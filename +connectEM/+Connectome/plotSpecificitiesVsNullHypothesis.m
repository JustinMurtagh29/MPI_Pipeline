% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'ProximalDendrite', 'PD'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'; ...
    'OtherDendrite', 'OD'};

targetLabels = transpose(targetClasses(:, 2));
targetClasses = transpose(targetClasses(:, 1));

% NOTE(amotta): Rows correspond to shaft and spine synapses, respectively.
% Columns correspond to the target classes listed above. All-zeros mark
% rates that were not measured (because biologically implausible).
inhFpRates = [ ...
     (1 / 20), (1 / 19), (5 / 17), (1 / 17), (4 / 18), (3 / 17);  % shaft
     (000000), (6 / 18), (6 / 18), (4 / 16), (000000), (6 / 17)]; % spine

minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);
[axonClasses.tag] = deal('Exc', 'Inh', 'TC', 'CC', 'Unclear');

%% build class connectome for shaft synapses
% This is needed to simulate the effect of FP inhibitory synapses.
clear cur*;

curShaftConn = conn;
curShaftConn.connectome.synIdx = cellfun( ...
    @(ids) ids(syn.synapses.type(ids) == 'Shaft'), ...
    curShaftConn.connectome.synIdx, 'UniformOutput', false);

curShaftConn = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        curShaftConn, axonClasses, 'minSynPre', minSynPre);
classConnShafts = ...
    connectEM.Connectome.buildClassConnectome( ...
        curShaftConn, 'targetClasses', targetClasses);
    
%% build complete class connectome
clear cur*;

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

classConn = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);
    
%% prepare for analysis while accounting for false positive inh. synapses

[axonClasses.synFalsePosRates] = deal([]);
[axonClasses.synFalsePosMethod] = deal([]);

axonClasses(2).synFalsePosRates = inhFpRates;
axonClasses(2).synFalsePosMethod = 'binomial';

%% plot
clear cur*;

evalT = cell(size(axonClasses));
barFracs = cell(size(axonClasses));

for curClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curClassId);
    
    curFpRates = curAxonClass.synFalsePosRates;
    curFpMethod = curAxonClass.synFalsePosMethod;
    
    % NOTE(amotta): If empirically determined synapse false positive rates
    % have been specified, we simulate the effect of synapse removal under
    % a binomial distribution. This necessitates multiple multiple runs to
    % get a feeling for the variability across trials.
    switch isempty(curFpRates)
        case true, numRuns = 1;
        case false, numRuns = 20;
    end
    
    curEvalT = array2table( ...
        nan(numRuns, 1 + numel(targetClasses)), ...
        'VariableNames', [{'Overall'}, targetClasses]);
    curBarFracs = nan(numRuns, 2 * numel(targetClasses) - 1);
    
    for curRun = 1:numRuns
        curFdrThresh = 0.2;
        curConn = classConn;
        
        if ~isempty(curFpRates)
            % Take into account false positive synapse rates
            curShaftConn = classConnShafts;
            curSpineConn = curConn - curShaftConn;

            rng(curRun);
            
            switch curFpMethod
              case 'binomial'
                % Binomial correction for false positive synapses
                curShaftConn = curShaftConn - binornd(curShaftConn, ...
                    repelem(curFpRates(1, :), size(curShaftConn, 1), 1));
                curSpineConn = curSpineConn - binornd(curSpineConn, ...
                    repelem(curFpRates(2, :), size(curSpineConn, 1), 1));
            
              case 'constFraction'
                % NOTE(amotta): Correction for false positive synapses
                % by removal of a random subset of synapses. With the
                % binomial approach the number of remove synapses
                % changes across runs. This method here keeps the
                % fraction of removed synapses constant.
                curShaftConn = ...
                    removeConstantFractionOfSynapses( ...
                        curShaftConn, curFpRates(1, :));
                curSpineConn = ...
                    removeConstantFractionOfSynapses( ...
                        curSpineConn, curFpRates(2, :));

              otherwise
                error('Invalid method "%s"', curFpMethod);
            end
            
            curConn = curShaftConn + curSpineConn;
        end
        
        curAxonClass.nullTargetClassProbs = ...
            connectEM.Specificity.calcFirstHitProbs( ...
                curConn(curAxonClass.nullAxonIds, :), ...
                'oneVersusRestBinomial');
        curAxonClass.specs =  ...
            connectEM.Specificity.calcForAxonClass( ...
                curConn, targetClasses, curAxonClass, ...
                'fdrThresh', curFdrThresh, ...
                'showPlot', false, ...
                'info', info);
       [curA, curB, curC] = ...
            connectEM.Specificity.calcAxonFractions( ...
                curAxonClass, targetClasses);
        
        curEvalT{curRun, :} = [curA, curB];
        curBarFracs(curRun, :) = curC;
    end
    
    evalT{curClassId} = curEvalT;
    barFracs{curClassId} = curBarFracs;
    
    curStatT = curEvalT;
    curStatT = varfun(@(v) [prctile(v, [0; 50; 100]); mean(v)], curStatT);
    curStatT.Properties.VariableNames = curEvalT.Properties.VariableNames;
    curStatT.Properties.RowNames = {'Min', 'Median', 'Max', 'Mean'};
    
    fprintf('%s\n\n', curAxonClass.title);
    disp(curStatT); fprintf('\n');
end

%% quantitative evaluation
clear cur*;

for curIdx = 1:numel(axonClasses)
    curAxonClass = axonClasses(curIdx);
    
    curEvalT = evalT{curIdx};
    curEvalT.fdrThresh = fdrThreshs(:);
    curEvalT = circshift(curEvalT, 1, 2);
    
    fprintf('%s\n\n', curAxonClass.title);
    disp(curEvalT); fprintf('\n');
end

%% Plot fraction of axons specific
% Copy & paste from +connectEM/+Figure/fractionOfAxonsSpecific.m
clear cur*;
curPlotClasses = 1:2;

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

for curIdx = 1:numel(curPlotClasses)
    curId = curPlotClasses(curIdx);
    curFracs = evalT{curId}.Overall;
    
    plot(curAx, repelem(curIdx, numel(curFracs)), curFracs');
end

curLines = curAx.Children;
set(curLines, 'Color', 'black', 'Marker', '.', 'MarkerSize', 16);

curAx.XLim = [0.5, numel(curPlotClasses) + 0.5];
curAx.YLim = [0, 1];

curAx.XTick = 1:numel(curPlotClasses);
curAx.XTickLabel = {axonClasses(curPlotClasses).tag};
ylabel(curAx, {'Fraction of'; 'axons specific'});

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [125, 160];

%% Plot distribution of specificities over target classes
clear cur*;

curBarWidth = 0.8;
curHistWidth = 0.09;
curPlotClasses = 1:2;

curColors = get(groot, 'defaultAxesColorOrder');
curColors = flip(curColors(1:numel(targetClasses), :), 1);

curMixed = (curColors(1:(end - 1), :) .^ 2 + curColors(2:end, :) .^ 2) / 2;
curMixed = cat(1, sqrt(curMixed), zeros(1, 3));

curColors = cat(1, transpose(curColors), transpose(curMixed));
curColors = transpose(reshape(curColors, 3, []));
curColors = num2cell(curColors(1:(end - 1), :), 2);

for curIdx = 1:numel(fdrThreshs)
    curFdrThresh = fdrThreshs(curIdx);
    
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [280, 280];

    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    curAx.YAxisLocation = 'right';

    curPlotData = cell2mat(cellfun( ...
        @(b) b(curIdx, :), barFracs(:), 'UniformOutput', false));
    curPlotData = flip(curPlotData(curPlotClasses, :), 2);

    curBars = bar(curAx, curPlotData, 'stacked', 'BarWidth', curBarWidth);
    set(curBars, {'FaceColor'}, curColors, 'EdgeColor', 'none');

    xlim(curAx, [0.5, numel(curPlotClasses) + 0.5]);
    xticklabels(curAx, {axonClasses(curPlotClasses).tag});
    xticks(curAx, 1:numel(curPlotClasses));

    curLeg = legend(curAx, ...
        flip(curBars(1:2:end)), targetLabels, ...
        'Location', 'EastOutside');
    title(curAx, sprintf('FDR threshold of %g', curFdrThresh));

    connectEM.Figure.config(curFig, info);
end

%% plotting
function classConn = ...
        removeConstantFractionOfSynapses(classConn, removeRate)
    assert(all(0 <= removeRate(:) & removeRate(:) <= 1));
    classCount = size(classConn, 2);
    
    if isscalar(removeRate)
        removeRate = repelem(removeRate, 1, classCount);
    else
        % NOTE(amotta): This will (intentionally) throw an error if the
        % number of elements in `removeRate` does not agree with the number
        % of target classes.
        removeRate = reshape(removeRate, 1, classCount);
    end
    
    for curIdx = 1:classCount
        curRemoveRate = removeRate(curIdx);
        
        % Find non-empty entries in class connectome
       [curIds, ~, curCounts] = find(classConn(:, curIdx));
        curIds = repelem(curIds(:), curCounts(:));
        
        % Select subset of synapses to drop
        curRemoveCount = round(curRemoveRate * numel(curIds));
        curIds = curIds(randperm(numel(curIds), curRemoveCount));
        
        % Update class connectome
       [curIds, ~, curCounts] = unique(curIds);
        curCounts = accumarray(curCounts, 1);
        
        classConn(curIds, curIdx) = ...
            classConn(curIds, curIdx) - curCounts;
    end
end
