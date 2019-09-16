% This script is based on 
%   +connectEM/+Connectome/plotGeometricPredictability.m
%   commit 8b9c96a8b57181c33f73cd411cf908138a76ca89
%
%   +connectEM/+Figure/coinnervation.m
%   commit f2418a0f552da11ce154faeffd01d65a4d867a11
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

[~, lengthFile] = fileparts(connFile);
lengthFile = sprintf('%s_pathLengths.mat', lengthFile);
lengthFile = fullfile(fileparts(connFile), lengthFile);

minSynPre = 10;
synTypes = {'All'};

% NOTE(amotta): Precision and recall values of automated synapse detection
% for correction of synapse frequencies / densities.
%   The precision and recall values for spine and shaft synapses were
% extracted from Supplementary Table 2. The overall precision / recall is
% used for spine synapses.
%   To disable this correction, unset `synPrecRecs`.
% synPrecRecs = [ ...
%     0.94, 0.89;  % spine
%     0.92, 0.69]; % shaft

axonClasses = { ...
    'Thalamocortical', 'TC'; ...
    'Corticocortical', 'CC'; ...
    'Inhibitory', 'Inh'};

targetClasses = { ...
    'ProximalDendrite', 'PD'; ...
    'SmoothDendrite', 'SD'; ...
    'ApicalDendrite', 'AD'; ...
    'OtherDendrite', 'Other'};

axonTags = axonClasses(:, 2);
axonClasses = axonClasses(:, 1);

targetTags = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

% Utility
renorm = @(v) v / sum(v);

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = ...
    connectEM.Connectome.load(param, connFile);
conn = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, [], 'minSynPre', minSynPre);

conn.denMeta.synCount = accumarray( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx), ...
   [numel(conn.dendrites), 1], [], 0);

lengths = load(lengthFile);

trunks = load(lengths.info.param.trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

maxSegId = Seg.Global.getMaxSegId(param);

%% Find trunk length for dendrite in connectome
clear cur*;

curLUT = Agglo.buildLUT(maxSegId, trunks);
conn.denMeta.trunkId = cellfun(@(ids) ...
    mode(nonzeros(curLUT(ids))), conn.dendrites);

curLens = lengths.trunkPathLengths;
curMask = ~isnan(conn.denMeta.trunkId);
conn.denMeta.trunkLength = nan(height(conn.denMeta), 1);
conn.denMeta.trunkLength(curMask) = curLens(conn.denMeta.trunkId(curMask));

%% Ignore axons and dendrites with too few synapses
% TODO(amotta): Should we also ignore postsynaptic targets?
clear cur*;

curAxonClasses = conn.axonMeta.axonClass;
curLastCat = categories(curAxonClasses);
curLastCat = curLastCat{end};

% NOTE(amotta): Add a new category called "Ignore" to the end of the
% list of all categories and mark all neurites with less than `minSyn`
% synapses as such.
curAxonClasses = addcats(curAxonClasses, {'Ignore'}, 'After', curLastCat);
curAxonClasses(conn.axonMeta.synCount < minSynPre) = 'Ignore';
conn.axonMeta.axonClass = curAxonClasses;

withoutIgnore = @(c) setdiff(c, {'Ignore'}, 'stable');
allAxonClasses = withoutIgnore(categories(conn.axonMeta.axonClass));
allTargetClasses = withoutIgnore(categories(conn.denMeta.targetClass));

%% Prepare synapse fractions
clear cur*;

curClassConn = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'axonClasses', allAxonClasses, ...
        'targetClasses', allTargetClasses);

[~, curAxonClassIds] = ismember(axonClasses, allAxonClasses);
[~, curTargetClassIds] = ismember(targetClasses, allTargetClasses);
curClassConn = curClassConn(curAxonClassIds, curTargetClassIds);

preSynSynFracs = reshape(renorm(sum(curClassConn, 2)), [], 1);
postSynSynFracs = reshape(renorm(sum(curClassConn, 1)), [], 1);

%% Prepare geometric data
% NOTE(amotta): The path lengths of soma-seeded agglomerates reported in
% the paper were derived from the ground truth skeleton tracings. Here,
% we're using the automatically calculated values.
clear cur*;

curLens = lengths.axonPathLengths;
curLens(conn.axonMeta.axonClass == 'Ignore') = 0;
[curMask, curIds] = ismember(conn.axonMeta.axonClass, axonClasses);
curFracs = accumarray(curIds(curMask), curLens(curMask));
preSynLengthFracs = curFracs / sum(curFracs);

curLens = conn.denMeta.trunkLength;
curLens(ismember(conn.denMeta.targetClass, {'Somata', 'Ignore'})) = 0;
[curMask, curIds] = ismember(conn.denMeta.targetClass, targetClasses);
curLens = accumarray(curIds(curMask), curLens(curMask));
postSynLengthFracs = curLens / sum(curLens);

%% Run analysis
clear cur*;

curPrecRecCorr = exist('synPrecRecs', 'var') && ~isempty(synPrecRecs);
if curPrecRecCorr; assert(isequal(size(synPrecRecs), [2, 2])); end
curPrecRecCorrName = {''; '; precision / recall corrected'};
curPrecRecCorrName = curPrecRecCorrName{1 + curPrecRecCorr};

curLogFac = @(ns) ...
    arrayfun(@(n) sum(log(1:n)), ns);
curMnLogLik = @(ns, ps) ...
    curLogFac(sum(ns, 1)) ...
  - sum(arrayfun(curLogFac, ns), 1) ...
  + sum(ns .* log(ps), 1);

for curSynIdx = 1%:numel(synTypes)
    curSynType = synTypes{curSynIdx};
    curSynLut = syn.synapses.type == curSynType;
    if strcmpi(curSynType, 'all'); curSynLut(:) = true; end
    
    curSpineLut = categories(syn.synapses.type);
    curSpineLut = curSpineLut(endsWith(curSpineLut, 'Spine'));
    curSpineLut = ismember(syn.synapses.type, curSpineLut);
    
    curConn = conn;
    curConn.connectome.synIdx = cellfun( ...
        @(ids) ids(curSynLut(ids)), ...
        curConn.connectome.synIdx, ...
        'UniformOutput', false);
    
    % NOTE(amotta): Correct for distortions in spine and shaft synapse
    % densities due to suboptimal precision / recall of automated synapse
    % detection.
    curSpineConn = curConn;
    curShaftConn = curConn;
   [curSpineConn.connectome.synIdx, ...
    curShaftConn.connectome.synIdx] = ...
        cellfun( ...
            @(ids) deal( ...
                ids( curSpineLut(ids)), ...
                ids(~curSpineLut(ids))), ...
            curConn.connectome.synIdx, ...
        'UniformOutput', false);
    
    curClassConns = cellfun(@(curConn) ...
        connectEM.Connectome.buildClassConnectome( ...
            curConn, ...
            'axonClasses', allAxonClasses, ...
            'targetClasses', allTargetClasses), ...
        {curSpineConn; curShaftConn}, ...
        'UniformOutput', false);
    
    if curPrecRecCorr
        curClassConns = cellfun( ...
            @(conn, precRec) conn * precRec(1) / precRec(2), ...
            curClassConns, num2cell(synPrecRecs, 2), ...
            'UniformOutput', false);
    end
    
    curClassConn = cat(3, curClassConns{:});
    curClassConn = sum(curClassConn, 3);
    curClassConn = round(curClassConn);
        
    curClassConnAbs = curClassConn;
    curSynCount = sum(curClassConn(:));
    curClassConn = curClassConn / curSynCount;
    
    curTitleStem = sprintf('%s (n = %d%s)', ...
        curSynType, curSynCount, curPrecRecCorrName);
    
    % Restrict to selected axon and target classes
   [~, curAxonClassIds] = ismember(axonClasses, allAxonClasses);
   [~, curTargetClassIds] = ismember(targetClasses, allTargetClasses);
    curClassConnAbs = curClassConnAbs(curAxonClassIds, curTargetClassIds);
    curClassConn = curClassConn(curAxonClassIds, curTargetClassIds);
    
    %% Statistical tests
    curMaxLogLik = @(pre, post) maxLogLikelihood( ...
        curClassConnAbs, preSynLengthFracs, postSynLengthFracs, pre, post);
   [curNullLogLik, curNullDof] = curMaxLogLik(false, false);
   
    fprintf('Statistical evaluation of Peters'' models\n');
    fprintf( ...
       ['* No synapse density correction: ', ...
        'log-likelihood of %g\n'], curNullLogLik);
       
    for curPreCorr = [false, true]
        for curPostCorr = [false, true]
            if ~curPreCorr && ~curPostCorr; continue; end
           [curLogLik, curDof] = curMaxLogLik(curPreCorr, curPostCorr);
           
            curChiDof = curDof - curNullDof;
            curLogLikRatio = -2 * (curNullLogLik - curLogLik);
            curPval = chi2cdf(curLogLikRatio, curChiDof, 'upper');
            
            curTitle = {'presynaptic', 'postsynaptic'};
            curTitle = curTitle([curPreCorr, curPostCorr]);
            curTitle = strjoin(curTitle, ' and ');
            curTitle(1) = upper(curTitle(1));
            
            fprintf( ...
               ['* %s density correction: ', ...
                'log-likelihood of %g, p-value of %g\n'], ...
                curTitle, curLogLik, curPval);
        end
    end
    
    %% Predicted synapse fraction is product of marginal synapse freqs.
    curTitle = strcat(curTitleStem, ' versus', ...
        ' product of synapse contributions');
    
    curAxonFracs = reshape(preSynSynFracs, [], 1);
    curTargetFracs = reshape(postSynSynFracs, 1, []);
    
    curRelClassConn = curClassConn ./ sum(curClassConn(:));
    curExpClassConn = curAxonFracs .* curTargetFracs;
    
    curLogLik = curMnLogLik(curClassConnAbs(:), curExpClassConn(:));
    fprintf('Log-likelihood for %s: %f\n', curTitle, curLogLik);
    
    curCorrCoeffs = curRelClassConn ./ curExpClassConn;
    curCorrCoeffs(curRelClassConn == 0 & curExpClassConn == 0) = 1;
    
    plotMatrix( ...
        info, curTitle, ...
        axonTags, curAxonFracs, ...
        targetTags, curTargetFracs, ...
        curRelClassConn, curExpClassConn, curCorrCoeffs);
    
    %% Predicted synapse fraction is product of marginal path length freqs.
    curTitle = strcat(curTitleStem, ' versus', ...
        ' product of path length contributions');
    
    curAxonFracs = reshape(preSynLengthFracs, [], 1);
    curTargetFracs = reshape(postSynLengthFracs, 1, []);
    
    curRelClassConn = curClassConn ./ sum(curClassConn(:));
    curExpClassConn = curAxonFracs .* curTargetFracs;
    
    curLogLik = curMnLogLik(curClassConnAbs(:), curExpClassConn(:));
    fprintf('Log-likelihood for %s: %f\n', curTitle, curLogLik);
    
    curCorrCoeffs = curRelClassConn ./ curExpClassConn;
    curCorrCoeffs(curRelClassConn == 0 & curExpClassConn == 0) = 1;
    
    plotMatrix( ...
        info, curTitle, ...
        axonTags, curAxonFracs, ...
        targetTags, curTargetFracs, ...
        curRelClassConn, curExpClassConn, curCorrCoeffs);
    
    %% Predicted synapse fraction is postsyn. marginal path length freqs.
    curTitle = strcat(curTitleStem, ' versus', ...
        ' postsynaptic path length contributions');
    
    curAxonFracs = reshape(preSynLengthFracs, [], 1);
    curTargetFracs = reshape(postSynLengthFracs, 1, []);
    
    curRelClassConn = curClassConn ./ sum(curClassConn, 2);
    curExpClassConn = repmat(curTargetFracs, size(curClassConn, 1), 1);
    
    curLogLik = sum(curMnLogLik(curClassConnAbs', curTargetFracs(:)));
    fprintf('Log-likelihood for %s: %f\n', curTitle, curLogLik);
    
    curCorrCoeffs = curRelClassConn ./ curExpClassConn;
    curCorrCoeffs(curRelClassConn == 0 & curExpClassConn == 0) = 1;
    
    plotMatrix( ...
        info, curTitle, ...
        axonTags, curAxonFracs, ...
        targetTags, curTargetFracs, ...
        curRelClassConn, curExpClassConn, curCorrCoeffs);
    
    %% Predicted synapse fraction is presyn. marginal path length freqs.
    curTitle = strcat(curTitleStem, ' versus', ...
        ' presynaptic path length contributions');
    
    curAxonFracs = reshape(preSynLengthFracs, [], 1);
    curTargetFracs = reshape(postSynLengthFracs, 1, []);
    
    curRelClassConn = curClassConn ./ sum(curClassConn, 1);
    curExpClassConn = repmat(curAxonFracs, 1, size(curClassConn, 2));
    
    curLogLik = sum(curMnLogLik(curClassConnAbs, curAxonFracs(:)));
    fprintf('Log-likelihood for %s: %f\n', curTitle, curLogLik);
    
    curCorrCoeffs = curRelClassConn ./ curExpClassConn;
    curCorrCoeffs(curRelClassConn == 0 & curExpClassConn == 0) = 1;
    
    plotMatrix( ...
        info, curTitle, ...
        axonTags, curAxonFracs, ...
        targetTags, curTargetFracs, ...
        curRelClassConn, curExpClassConn, curCorrCoeffs);
end

%% Statistical tests
function [maxLogLik, numDof] = maxLogLikelihood( ...
        classConn, preNullFracs, postNullFracs, corrPre, corrPost)
    renorm = @(v) v / sum(v(:));
    
    preN = size(classConn, 1);
    postN = size(classConn, 2);
    
    preVars = corrPre * preN;
    postVars = corrPost * postN;
    
    % NOTE(amotta): Check if we can omit calculation of the first two terms
    % in `curMnLogLik`. They seem to be independent of the probabilities.
    logFac = @(ns) ...
        arrayfun(@(n) sum(log(1:n)), ns);
    mnLogLik = @(ns, ps) ...
        logFac(sum(ns, 1)) ...
      - sum(arrayfun(logFac, ns), 1) ...
      + sum(ns .* log(ps), 1);
    
    ps = @(x) renorm(reshape( ...
        (preNullFracs(:) .* [ ...
            x(1:preVars); ones(preN - preVars, 1)]) ...
     .* (postNullFracs(:) .* [ ...
            x((preVars + 1):end); ...
            ones(postN - postVars, 1)])', ...
        [], 1));
    
    if corrPre || corrPost
        x0 = [ ...
            renorm(ones(preVars, 1)); ...
            renorm(ones(postVars, 1))];

        Aeq = [ ...
            ones(corrPre, preVars), zeros(corrPre, postVars);
            zeros(corrPost, preVars), ones(corrPost, postVars)];
        beq = ones(corrPre + corrPost, 1);

        lb = zeros(preVars + postVars, 1);
        ub = ones(preVars + postVars, 1);
        
        opts = optimoptions('fmincon', 'Display', 'notify');

       [~, minNegLogLik] = fmincon( ...
            @(x) -mnLogLik(classConn(:), ps(x)), ...
            x0, [], [], Aeq, beq, lb, ub, [], opts);
        maxLogLik = -minNegLogLik;
    else
        maxLogLik = mnLogLik(classConn(:), ps([]));
    end
    
    numDof = ...
        corrPre * (preN - 1) ...
      + corrPost * (postN - 1);
end

%% Plotting
function plotMatrix( ...
        info, titleStr, ...
        axonTags, axonFracs, ...
        targetTags, targetFracs, ...
        relClassConn, expClassConn, corrCoeffs)
    % Configuration
    minObsExpRatio = log10(1.25);
    colorLim = log10(5);
    
    rows = numel(axonTags);
    cols = numel(targetTags);
    frac = rows / cols;
    
    colorMat = log10(corrCoeffs);
    colorMap = connectEM.Figure.redBlue(129);
    
    curVals = linspace(-colorLim, +colorLim, size(colorMap, 1));
    colorMap(abs(curVals) < minObsExpRatio, :) = 1;

    fig = figure();
    ax = axes(fig);

    imshow( ...
        colorMat, [-colorLim, +colorLim], ...
        'Colormap', colorMap, ...
        'Parent', ax);

    fig.Color = 'white';
    fig.Position(3:4) = 1000 .* [1, frac];

    ax.Visible = 'on';
    ax.TickDir = 'out';
    ax.Box = 'off';

    ax.XAxisLocation = 'top';
    ax.XTick = 1:size(colorMat, 2);
    ax.XTickLabel = arrayfun( ...
        @(t, p) sprintf('%s (%.2g)', t{1}, 100 * p), ...
        targetTags(:), targetFracs(:), 'UniformOutput', false);
    ax.XTickLabelRotation = 90;

    ax.YTick = 1:size(colorMat, 1);
    ax.YTickLabel = arrayfun( ...
        @(t, p) sprintf('%s (%.2g)', t{1}, 100 * p), ...
        axonTags(:), axonFracs(:), 'UniformOutput', false);
    ax.Position = [0.1, 0.01, 0.75, 0.75];

    for curIdx = 1:numel(colorMat)
       [curRow, curCol] = ind2sub(size(colorMat), curIdx);
        curEdgeColor = 'none';

        curBoxSize = ax.Position(3:4) ./ [cols, rows];
        curOff = [curCol, numel(axonTags) - curRow + 1];
        curOff = ax.Position(1:2) + (curOff - 1) .* curBoxSize;

        curAnn = annotation(fig, ...
            'textbox', [curOff, curBoxSize], ...
            'String', { ...
                ratio2str(corrCoeffs(curIdx), 3); ...
                sprintf('(%.2g vs. %.2g)', ...
                    100 * relClassConn(curIdx), ...
                    100 * expClassConn(curIdx))});

        curAnn.HorizontalAlignment = 'center';
        curAnn.VerticalAlignment = 'middle';
        curAnn.EdgeColor = curEdgeColor;
        curAnn.Color = 'black';
        curAnn.FontSize = 12;
        curAnn.LineWidth = 2;
    end

    cbar = colorbar('peer', ax);
    cbar.Label.String = '"Affinity" (measured / expected)';
    cbar.Ticks = [0, minObsExpRatio, colorLim];
    cbar.Ticks = unique(cat(2, -cbar.Ticks, +cbar.Ticks));
    cbar.TickLabels = arrayfun( ...
        @(f) ratio2str(10 ^ f), ...
        cbar.Ticks, 'UniformOutput', false);
    cbar.TickDirection = 'out';
    cbar.Position = [0.86, 0.1, 0.02, 0.8];
    cbar.Position([2, 4]) = ax.Position([2, 4]);

    title(ax, ...
        {info.filename; info.git_repos{1}.hash; titleStr}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

function str = ratio2str(frac, prec)
    fmt = '%s%g';
    if exist('prec', 'var') && ~isempty(prec)
        fmt = ['%s%.', num2str(prec, '%d'), 'g'];
    end
    
    str = {'/', '', '×'};
    str = str{2 + sign(log(frac))};
    if frac < 1; frac = 1 / frac; end
    str = sprintf(fmt, str, frac);
end
