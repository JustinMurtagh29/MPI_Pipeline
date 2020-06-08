% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connNames = { ...
    'InitialAxonAgglomerates', ...
    'FinalAxonAgglomerates'};
connFiles = { ...
    'connectome_axons-04_dendrites-wholeCells-autoSpines-v1-classified-v2_SynapseAgglos-autoPreRobo-v1-classified.mat';
    'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat'};
connFiles = fullfile(rootDir, 'connectomeState', connFiles);

% NOTE(amotta): For evaluation of spine heads
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
clear cur*;

param = fullfile(rootDir, 'allParameter.mat');
param = Util.load(param, 'p');

maxSegId = Seg.Global.getMaxSegId(param);
curShAgglos = Util.load(shFile, 'shAgglos');

[conns, syns] = deal(cell(size(connFiles)));
for curConnIdx = 1:numel(connFiles)
    curConnFile = connFiles{curConnIdx};
    curConn = load(curConnFile);
    
    curSynFile = curConn.info.param.synFile;
    curSyn = load(curSynFile);
    
    % Complete axon meta data (for primary spine synapse fraction)
    curConn.axonMeta = completeAxonMeta(param, curConn, curSyn);
    
    % Complete dendrite meta data (number of spine heads)
    curDendLUT = Agglo.buildLUT(maxSegId, curConn.dendrites);
    
    curShDendIds = cellfun( ...
        @(segIds) mode(nonzeros(curDendLUT(segIds))), ...
        curShAgglos, 'UniformOutput', false);
    
    assert(all(cellfun(@isscalar, curShDendIds)));
    curShDendIds = cell2mat(curShDendIds(:));
    
    % NOTE(amotta): `mode` returns `NaN` on empty inputs. For our count,
    % let's remove spine heads that did not get attached.
    curShDendIds(isnan(curShDendIds)) = [];
    assert(all(curShDendIds(:) > 0));
    
    curConn.denMeta.spineHeadCount = accumarray( ...
        curShDendIds(:), 1, size(curConn.dendrites(:)));
    
    conns{curConnIdx} = curConn;
    syns{curConnIdx} = curSyn;
end

%% collect numbers
clear cur*;
connVals = cell(numel(conns), 2);

for idx = 1:numel(conns)
    vals = cell(0, 2);
    conn = conns{idx};
    
   [~, connName] = fileparts(connFiles{idx});
    vals(end + 1, :) = {'Connectome', connName}; %#ok
    
    vals(end + 1, :) = { ...
        '# axons', numel(conn.axons)}; %#ok
    vals(end + 1, :) = { ...
        '# postsynaptic targets', numel(conn.dendrites)}; %#ok
    vals(end + 1, :) = { ...
        '# spine heads', sum(conn.denMeta.spineHeadCount)}; %#ok
    
    vals(end + 1, :) = { ...
        '# synapses', conn.connectomeMeta.noSynapses}; %#ok
    vals(end + 1, :) = { ...
        '# spine synapses', sum(conn.axonMeta.spineSynCount)}; %#ok
    
    vals(end + 1, :) = { ...
        '# axons without synapses', ...
        sum(conn.axonMeta.synCount == 0)}; %#ok
    vals(end + 1, :) = { ...
        '# axons with at least 10 synapses', ...
        sum(conn.axonMeta.synCount >= 10)}; %#ok

    postSynCount = accumarray( ...
        conn.connectome.edges(:, 2), ...
        cellfun(@numel, conn.connectome.synIdx)', ...
       [numel(conn.dendrites), 1]);
    
    vals(end + 1, :) = { ...
        '# targets without synapses', ...
        sum(postSynCount == 0)}; %#ok
    vals(end + 1, :) = { ...
        '# targets with at least 10 synapses', ...
        sum(postSynCount >= 10)}; %#ok
    
    for classIdx = 1:numel(conn.denClasses)
        vals(end + 1, :) = { ...
            sprintf('# synapses onto %s', conn.denClasses{classIdx}), ...
            sum(conn.classConnectome(:, classIdx))}; %#ok
    end
    
    connVals{idx, 1} = vals(:, 1);
    connVals{idx, 2} = vals(:, 2);
end

%% Show results
rowNames = cat(1, connVals{:, 1});
rowNames = unique(rowNames, 'stable');

tableData = cell(numel(rowNames), numel(conns));
tableData(:) = num2cell(nan);

for idx = 1:numel(conns)
   [~, rows] = ismember(connVals{idx, 1}, rowNames);
    tableData(rows, idx) = connVals{idx, 2};
end

tableData(1, :) = [];
rowNames(1, :) = [];

t = cell2table(tableData);
t.Properties.VariableNames = connNames;
t.Properties.RowNames = rowNames;

format('long', 'g');
disp(t);

%% Plot primary spine synapse fraction distribution
% See +connectEM/+Figure/spineSynapseFraction.m
clear cur*;
curMinSynCount = 10;
curThreshLines = [0.2, 0.5];
curBinEdges = linspace(0, 1, 21);

curFig = figure();
curAx = axes(curFig);
axis(curAx, 'square')
hold(curAx, 'on');

for curConnIdx = 1:numel(connFiles)
    curConn = conns{curConnIdx};
    curAxonMeta = curConn.axonMeta;

    % Remove axons with too few synapses
    curAxonMeta(curAxonMeta.synCount < curMinSynCount, :) = [];

    curAxonMeta.fullPriSpineSynFrac = ...
        curAxonMeta.fullPriSpineSynCount ...
     ./ curAxonMeta.fullSynCount;
 
    histogram(curAx, ...
        curAxonMeta.fullPriSpineSynFrac, curBinEdges, ...
        'Normalization', 'probability');
end

lines = arrayfun(@(x) plot(curAx, [x, x], curAx.YLim), curThreshLines);
set(lines, 'Color', 'black', 'LineWidth', 2', 'LineStyle', '--');

xlabel(curAx, {'Fraction of synapses'; 'onto primary spines'});
ylabel(curAx, 'Fraction of axons');

curLeg = legend(curAx, connNames);
curLeg.Location = 'EastOutside';

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [500, 265];

%% Utilities
function axonMeta = completeAxonMeta(param, conn, syn)
    % NOTE(amotta): This funtion was extracted form
    % +connectEM/+Axon/completeSynapseMeta.m
    synapses = syn.synapses;
    synapses.id = reshape( ...
        1:size(synapses, 1), [], 1);
    synapses.ontoSpine = ...
        synapses.type == 'PrimarySpine' ...
      | synapses.type == 'SecondarySpine';
  
    maxSegId = Seg.Global.getMaxSegId(param);
    axonLUT = Agglo.buildLUT(maxSegId, conn.axons);
    
    synapses.axonId = cellfun( ...
        @(segIds) setdiff(axonLUT(segIds), 0), ...
        synapses.presynId, 'UniformOutput', false);
    
    synapses(~cellfun(@isscalar, synapses.axonId), :) = [];
    synapses.axonId = cell2mat(synapses.axonId);

    axonMeta = conn.axonMeta;
    axonMeta.fullSynCount = accumarray( ...
        synapses.axonId, 1, size(axonMeta.id));
    axonMeta.fullSpineSynCount = accumarray( ...
        synapses.axonId, synapses.ontoSpine, size(axonMeta.id));
    priSynSpineMask = ...
        synapses.type == 'PrimarySpine';
    axonMeta.fullPriSpineSynCount = accumarray( ...
        synapses.axonId, priSynSpineMask, size(axonMeta.id));
end
