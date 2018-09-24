% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSynPre = 10;
minSynPost = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

[~, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

%% Primary spine synapse fraction
clear cur*;

[curTypes, ~, curTemp] = unique(syn.synapses.type);
curTemp = accumarray(curTemp, 1);

% Fraction of all spine being primary spine innervations
priSpineSynFrac = curTemp(curTypes == 'PrimarySpine') / sum(curTemp) %#ok

%% Synapse type versus polarity
clear cur*;
    
curLUT = zeros(maxSegId, 1);
curLUT(cell2mat(conn.axons(axonClasses(1).axonIds))) = +1;
curLUT(cell2mat(conn.axons(axonClasses(2).axonIds))) = -1;

syn.synapses.polarities = cellfun( ...
    @(ids) reshape(setdiff(curLUT(ids), 0), [], 1), ...
    syn.synapses.presynId, 'UniformOutput', false);

curT = table;
curT.type = repelem(syn.synapses.type, ...
    cellfun(@numel, syn.synapses.polarities));
curT.polarity = cell2mat(syn.synapses.polarities);

synTypes = table;
synTypes.type = unique(curT.type);
synTypes.excInhRatio = arrayfun(@(t) ...
    sum(curT.type == t & curT.polarity > 0) ...
  / sum(curT.type == t), synTypes.type);
disp(synTypes)

%% Excitatory / inhibitory fraction over volume
clear cur*;
    
curExcAxonIds = axonClasses(1).axonIds;
curExcSynCount = sum(conn.axonMeta.fullSynCount(curExcAxonIds));

curInhAxonIds = axonClasses(2).axonIds;
curInhSynCount = sum(conn.axonMeta.fullSynCount(curInhAxonIds));

volExcInhRatio = curExcSynCount / (curExcSynCount + curInhSynCount) %#ok

%% Synapse density for excitatory / inhibitory axons
clear cur*;

curDensity = conn.axonMeta.fullSynCount ./ conn.axonMeta.pathLen;
curAxonClasses = {axonClasses(1).axonIds; axonClasses(2).axonIds};

axonTypes = table;
axonTypes.type = categorical({'excitatory'; 'inhibitory'});
axonTypes.axonCount = cellfun(@numel, curAxonClasses);
axonTypes.synDensityMean = cellfun( ...
    @(ids) mean(curDensity(ids)), curAxonClasses);
axonTypes.synDensityStd = cellfun( ...
    @(ids) std(curDensity(ids)), curAxonClasses);

axonTypes %#ok

%% Synapse density along dendrites
clear cur*;

curLUT = Agglo.buildLUT(maxSegId, conn.dendrites);

syn.synapses.dendIds = cellfun( ...
    @(ids) reshape(setdiff(curLUT(ids), 0), [], 1), ...
    syn.synapses.postsynId, 'UniformOutput', false);

curT = table;
curT.synId = transpose(repelem( ...
    1:numel(syn.synapses.dendIds), ...
    cellfun(@numel, syn.synapses.dendIds)));
curT.polarity = syn.synapses.polarities(curT.synId);
curT.dendId = cell2mat(syn.synapses.dendIds);

curAccum = @(subs, vals) accumarray(subs, vals, [height(conn.denMeta), 1]);
conn.denMeta.fullSynCount = curAccum(curT.dendId, 1);

curPolarity = cell2mat(curT.polarity);
curT = curT(repelem(1:height(curT), cellfun(@numel, curT.polarity)), :);
curT.polarity = curPolarity;

conn.denMeta.excSynCount = curAccum(curT.dendId, curT.polarity > 0);
conn.denMeta.inhSynCount = curAccum(curT.dendId, curT.polarity < 0);
conn.denMeta.classSynCount = curAccum(curT.dendId, 1);

conn.denMeta.excInhRatio = ...
    conn.denMeta.excSynCount ./ ( ...
    conn.denMeta.excSynCount + conn.denMeta.inhSynCount);

meanDendExcInhRatio = mean( ...
    conn.denMeta.excInhRatio( ...
        conn.denMeta.classSynCount >= minSynPost)) %#ok

dendTypes = table;
[dendTypes.type, ~, curTemp] = unique(conn.denMeta.targetClass);

dendTypes.dendIds = accumarray( ...
    curTemp, reshape(1:height(conn.denMeta), [], 1), ...
    [], @(ids) {ids(conn.denMeta.classSynCount(ids) >= minSynPost)});
dendTypes.dendCount = cellfun(@numel, dendTypes.dendIds);

dendTypes.excInhRatioMean = cellfun( ...
    @(ids) mean(conn.denMeta.excInhRatio(ids)), dendTypes.dendIds);
dendTypes.excInhRatioStd = cellfun( ...
    @(ids) std(conn.denMeta.excInhRatio(ids)), dendTypes.dendIds);

dendTypes %#ok

%% Targets of excitatory shaft synapses
clear cur*;

curExcSyns = syn.synapses;
curExcSyns = curExcSyns(cellfun( ...
    @(p) isequal(p, 1), curExcSyns.polarities), :);

curExcShaftSyns = curExcSyns;
curExcShaftSyns = curExcShaftSyns( ...
    curExcShaftSyns.type == 'Shaft', :);

excShaftSyns = table;
excShaftSyns.targetClass = unique(conn.denMeta.targetClass);

curTemp = conn.denMeta.targetClass(cell2mat(curExcShaftSyns.dendIds));
[~, curTemp] = ismember(curTemp, excShaftSyns.targetClass);

excShaftSyns.excShaftFrac = accumarray(curTemp, 1);
excShaftSyns.excShaftFrac = ...
    excShaftSyns.excShaftFrac ...
  / sum(excShaftSyns.excShaftFrac);


curTemp = conn.denMeta.targetClass(cell2mat(curExcSyns.dendIds));
[~, curTemp] = ismember(curTemp, excShaftSyns.targetClass);

excShaftSyns.excFrac = accumarray(curTemp, 1);
excShaftSyns.excFrac = excShaftSyns.excFrac / sum(excShaftSyns.excFrac);

excShaftSyns %#ok
