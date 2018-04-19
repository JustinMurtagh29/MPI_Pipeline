% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

dendFile = fullfile(rootDir, 'aggloState', 'dendrites_18.mat');
wcFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4.mat');
outFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_v6.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

dend = load(dendFile);
wc = load(wcFile);

%% Build output
dendAgglos = dend.dendrites;
wcAgglos = wc.wholeCellsNoAxon(:);

out = struct;
out.dendrites = [dendAgglos; wcAgglos];
out.indBigDends = [dend.indBigDends; true(size(wcAgglos))];

% Whole cell indices
out.indWholeCells = [false(size(dendAgglos)); true(size(wcAgglos))];
out.idxWholeCells = [zeros(size(dendAgglos)); (1:numel(wcAgglos))'];

% Axon initial segment indices
out.indAIS = [dend.indAIS; false(size(wcAgglos))];
out.idxAIS = [dend.idxAIS; zeros(size(wcAgglos))];

out.info = info;

%% Mask overlapping segments
% Nodes with segments IDs that are both in the dendrite and whole cell
% agglomerates will be mapped to segment ID NaN.
maxSegId = Seg.Global.getMaxSegId(param);
maxSegId = maxSegId + 1;

maskLUT = ...
    Agglo.buildLUT(maxSegId, Agglo.fromSuperAgglo(wcAgglos)) ...
  & Agglo.buildLUT(maxSegId, Agglo.fromSuperAgglo(dendAgglos));
maskCount = sum(maskLUT);

if maskCount > 0
    warning('Masking away %d overlapping segments', maskCount);
    
    for curIdx = 1:numel(out.dendrites)
        curDend = out.dendrites(curIdx);
        curMask = curDend.nodes(:, 4);
        curMask(isnan(curMask)) = maxSegId;
        curMask = maskLUT(curMask);
        
        curDend.nodes(curMask, 4) = nan;
        out.dendrites(curIdx) = curDend;
    end
end

%% Sanity check
out.dendrites = SuperAgglo.clean(out.dendrites);

%% Writing result
Util.saveStruct(outFile, out);
Util.protect(outFile);
